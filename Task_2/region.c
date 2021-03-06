#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "bmp.h"

typedef struct{
    int x;
    int y;
} pixel_t;

typedef struct{
    int size;
    int buffer_size;
    pixel_t* pixels;
} stack_t;


// Global variables
int rank,                       // MPI rank
    size,                       // Number of MPI processes
    dims[2],                    // Dimensions of MPI grid
    coords[2],                  // Coordinate of this rank in MPI grid
    periods[2] = {0,0},         // Periodicity of grid
    north,south,east,west,      // Four neighbouring MPI ranks
    image_size[2] = {512,512},  // Hard coded image size
    local_image_size[2];        // Size of local part of image (not including border)

	
MPI_Comm cart_comm;             // Cartesian communicator


// MPI datatypes, you may have to add more.
MPI_Datatype border_row_t,
             border_col_t,
             original_image,
             partial_image,
			 partial_region;


unsigned char *image,           // Entire image, only on rank 0
              *region,          // Region bitmap. 1 if in region, 0 elsewise
              *local_image,     // Local part of image
              *local_region;    // Local part of region bitmap
			  
// Create new pixel stack
stack_t* new_stack(){
    stack_t* stack = (stack_t*)malloc(sizeof(stack_t));
    stack->size = 0;
    stack->buffer_size = 1024;
    stack->pixels = (pixel_t*)malloc(sizeof(pixel_t*)*1024);
}


// Push on pixel stack
void push(stack_t* stack, pixel_t p){
    if(stack->size == stack->buffer_size){
        stack->buffer_size *= 2;
        stack->pixels = realloc(stack->pixels, sizeof(pixel_t)*stack->buffer_size);
    }
    stack->pixels[stack->size] = p;
    stack->size += 1;
}


// Pop from pixel stack
pixel_t pop(stack_t* stack){
    stack->size -= 1;
    return stack->pixels[stack->size];
}
 

// Check if two pixels are similar. The hardcoded threshold can be changed.
// More advanced similarity checks could have been used.
int similar(unsigned char* im, pixel_t p, pixel_t q){
    int a = im[p.x +  p.y * image_size[1]];
    int b = im[q.x +  q.y * image_size[1]];
    int diff = abs(a-b);
    return diff < 2;
}


// Create and commit MPI datatypes
void create_types(){
    int count = 512 / dims[1];//this is true if the count in mpi_vector is defined for y-cord for our vector.
    int blocklen = 512/dims[0];// this is true if blocklen is for the block total lenght.
    int stride = 512 - blocklen + 1;
   
    MPI_Type_vector(count, blocklen, stride, MPI_INT, &partial_image );//here, the count would be the total amount of pixels divide by total amount of cores.
    MPI_Type_commit( &partial_image);
    MPI_Type_vector(count, 1, 512, MPI_INT, &border_col_t);
    MPI_Type_commit(&border_col_t);

    MPI_Type_vector(1, blocklen, 1, MPI_INT, &border_row_t );
    MPI_Type_commit(&border_row_t); 
	
}
// The initialization of the halo.
void init_halo(){
    int max_x_value = local_image_size[1] + 1; 
    int max_y_value = local_image_size[0] + 1;
    //each partial_image, will send their respective borders to their neighbour. This section is actually just an halo exchange
        if(north >= 0){
            MPI_Send(&local_image[local_image_size[1] + 1], 1, border_row_t,   north,  1,cart_comm);
            MPI_Recv(&local_image[0], 1, border_row_t,   north,  1, cart_comm, MPI_STATUS_IGNORE);
        }
        if(south >= 0){
            MPI_Send(&local_image[(max_y_value) * 1 + (local_image_size[0]-2)], 1, border_row_t,   south,  1, cart_comm);
            MPI_Recv(&local_image[max_x_value*local_image_size[0]-1], 1, border_row_t,   south,  1, cart_comm, MPI_STATUS_IGNORE);
        }
        if(east >= 0){
            MPI_Send(&local_image[local_image_size[1]-2], 1, border_col_t,   east,  1, cart_comm);
            MPI_Recv(&local_image[local_image_size[1]-1], 1, border_col_t,   east,  1, cart_comm, MPI_STATUS_IGNORE);
        }
        if(west >= 0){
            MPI_Send(&local_image[1], 1, border_col_t,   west,  1, cart_comm);
            MPI_Recv(&local_image[0], 1, border_col_t,   west,  1, cart_comm, MPI_STATUS_IGNORE);
        }

}

// Send image from rank 0 to all ranks, from image to local_image
// The function should react different acording to the rank of the program, rank 0 is the sender, rank 1-(size-1) is the recivers
void distribute_image(){
    // int local_image_x = 512/ dim[0];
    // int local_image_y = 512/ dim[1];
    // MPI_scatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator);  

    // My original idea was using MPI_scatter to distribute the partitioned image. But after discussing it with the T.A, I was informed that my idea was incorrect.

    //Here my thoughts are, I have allready created a MPI_vector in create_types. Then i can start using it here.
    if(rank ==0){
		int rank_counter=0;
        for(int row=0 ; row<dims[0]; row++){
            for(int col = 0; col < dims[1]; col ++ ){
			rank_counter+=1;
                MPI_Send(&image[(row * local_image_size[1]) + (col*local_image_size[0])], 1, partial_image, rank_counter, 1, cart_comm );
            }
       // MPI_Send(&local_image_x, 1, MPI_INT, i, 1);
       // MPI_Send(&local_image-y, 1, MPI_INT, i, 1);
        }
       //soem kind of algorithm to store the local element to rank 0, using the same logick as the MPI_vector, my initial idea was to send an mpi message to myself, but quickly gave up the idea.

	   
        for(int row = 0; row< local_image_size[1]; row++){
            for (int col; col < local_image_size[0]; col++){
                
                int index = (row)*local_image_size[1] + col;
                int destination = (row+1)*local_image_size[1] +col+1 ;
                local_image [destination] = image[index];
            }
        }

    }

    else{
        MPI_Recv(&local_image[local_image_size[0]*1 + 1], 1, partial_image, 0, 1, cart_comm, MPI_STATUS_IGNORE);
        // MPI_Recv(&local_image_size[0], 1, MPI_INT, 0, 1);
        // MPI_Recv(&local_image_size[1], 1, MPI_INT, 0, 1);

    }
    //Initiate the halo for all the process, doing so by asking every rank to send their border to their respectable naighbour.
    init_halo();
   // MPI_Scatter(image, 1 , partial_image, local_image, 1,partial_image, 0,cart_comm );
 
}


// Exchange borders with neighbour ranks
void exchange(stack_t* stack){  
    unsigned char *buffer =(unsigned char*) malloc(sizeof(unsigned char)*local_image_size[0] + 2); 
       //hvorfor trenger exchange stack?
        if(north){
            MPI_Send(&local_region[0], 1, border_row_t,   north, 1, cart_comm);
          
            MPI_Recv(&buffer, 1, border_row_t,   north, 1, cart_comm, MPI_STATUS_IGNORE);
            for(int i = 0; i< local_image_size[1]; i++ ){
                if(buffer[i]==1 && local_region[i] != 1){
                    pixel_t candidate;
                    candidate.x = i;
                    candidate.y = 0;
            
                    push(stack, candidate);
                }
            }
        }
        if(south){
            MPI_Send(&local_region[local_image_size[0] + 1], 1, border_row_t,   south,  1, cart_comm);
           
            MPI_Recv(&buffer, 1, border_row_t,   south,  1, cart_comm, MPI_STATUS_IGNORE);
            for(int i = 0; i<local_image_size[1]; i++){
                if(buffer[i] ==1 && local_region[ (local_image_size[1]*local_image_size[0])+i ] != 1){
                    pixel_t candidate;
                    candidate.x = i;
                    candidate.y = local_image_size[1];
            
                    push(stack, candidate);

                }

            }
        }
        if(east){
         
            MPI_Send(&local_region[local_image_size[1]+1] , 1, border_col_t,   east,  1, cart_comm);
          
            MPI_Recv(&buffer, 1, border_col_t,   east,  1, cart_comm, MPI_STATUS_IGNORE);
            for(int i = 0; i<local_image_size[1]; i++){
                pixel_t candidate;
                    candidate.x = local_image_size[0];
                    candidate.y = i;
            
                    push(stack, candidate);

            }

        }
        if(west){
            MPI_Send(&local_region[0], 1, border_col_t,   west,  1, cart_comm);
           
            MPI_Recv(&buffer, 1, border_col_t,   west,  1, cart_comm, MPI_STATUS_IGNORE);
            for(int i = 0; i<local_image_size[1]; i++){
                pixel_t candidate;
                    candidate.x = 0;
                    candidate.y = i;
            
                    push(stack, candidate);

            }
        } 
        free(buffer);
}


// Gather region bitmap from all ranks to rank 0, from local_region to region
void gather_region(){
    //Rank 0 will start collecting the others llocal image.
    if(0 == rank){
		int rank_counter =0;
		for(int row = 0; row< dims[1]; row++ ){
			for(int col = 1; col<dims[0]; col++){
				rank_counter += 1;
				MPI_Recv(&region[(row * local_image_size[1]) + (col*local_image_size[0])], 1, partial_image , rank_counter, 1, cart_comm, MPI_STATUS_IGNORE);
			}
		}
    }

    else{
        MPI_Send(&local_region, 1, partial_image, 0, 1,cart_comm);
        

    }
}

// Determine if all ranks are finished. You may have to add arguments.
// You dont have to have this check as a seperate function
// Since i put the exchange in the growe region, 
int finished(stack_t * stack){
   //the stack was emptied. Ergo this rank recived nothing from his/hers neighbour. The result is that this rank will broadcast his result to everyone. 
   // MPI_Bcast(&buf, 1, MPI_INT, root, MPI_COMM_WORLD)';
   int* finished = (int*)malloc(sizeof(int)*size);
   if(stack->size ==0){
		finished[rank] = 1;
   }
   else{
		finished[rank] = 0;
   }
   for(int i =0; i< size; i++){
		MPI_Bcast(&finished[i], 1, MPI_INT, finished[i], cart_comm);
   }
   for(int i =0; i<size; i++){
	if(finished[i]==0){
		return 0;
	}
	}
	return 1;
}


// Check if pixel is inside local image
int inside(pixel_t p){
    return (p.x >= 0 && p.x < local_image_size[1] && p.y >= 0 && p.y < local_image_size[0]);
}

// Adding seeds in corners.
void add_seeds(stack_t* stack){
    int seeds [8];
    seeds[0] = 5;
    seeds[1] = 5;
    seeds[2] = image_size[1]-5;
    seeds[3] = 5;
    seeds[4] = image_size[1]-5;
    seeds[5] = image_size[0]-5;
    seeds[6] = 5;
    seeds[7] = image_size[0]-5;
    
    for(int i = 0; i < 4; i++){
        pixel_t seed;
        seed.x = seeds[i*2] - coords[1]*local_image_size[1];
        seed.y = seeds[i*2+1] -coords[0]*local_image_size[0];
        
        if(inside(seed)){
            push(stack, seed);
        }
    }
}


// Region growing, serial implementation
void grow_region(){
    
    stack_t* stack = new_stack();
    add_seeds(stack);
    
    while(stack->size > 0){
        pixel_t pixel = pop(stack);
        
        region[pixel.y * image_size[1] + pixel.x] = 1;
        
        
        int dx[4] = {0,0,1,-1}, dy[4] = {1,-1,0,0};
        for(int c = 0; c < 4; c++){
            pixel_t candidate;
            candidate.x = pixel.x + dx[c];
            candidate.y = pixel.y + dy[c];
            
            if(!inside(candidate)){
                continue;
            }
            
            
            if(region[candidate.y * image_size[1] + candidate.x]){
                continue;
            }
            
            if(similar(image, pixel, candidate)){
                region[candidate.x + candidate.y * image_size[1]] = 1;
                push(stack,candidate);
            }
        }
        if(stack->size == 0){
            exchange(stack); // this line should be able to let the loop continue, If the stack is empty, it will call a exchange, after the exchange the stack-> size should be non-zero
                            // the next line when the while loop test it, it will show as true and therefor continue, with this, we would not need to know when 
        }

    }

}


// MPI initialization, setting up cartesian communicator
void init_mpi(int argc, char** argv){
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm );
    MPI_Cart_coords( cart_comm, rank, 2, coords );
    
    MPI_Cart_shift( cart_comm, 0, 1, &north, &south );
    MPI_Cart_shift( cart_comm, 1, 1, &west, &east );
}


void load_and_allocate_images(int argc, char** argv){

    if(argc != 2){
        printf("Useage: region file");
        exit(-1);
    }
    
    if(rank == 0){
        image = read_bmp(argv[1]);
        region = (unsigned char*)calloc(sizeof(unsigned char),image_size[0]*image_size[1]);
    }
    
    local_image_size[0] = image_size[0]/dims[0];
    local_image_size[1] = image_size[1]/dims[1];
    
    int lsize = local_image_size[0]*local_image_size[1];
    int lsize_border = (local_image_size[0] + 2)*(local_image_size[1] + 2);
    local_image = (unsigned char*)malloc(sizeof(unsigned char)*lsize_border);
    local_region = (unsigned char*)calloc(sizeof(unsigned char),lsize_border);
}


void write_image(){
    if(rank==0){
        for(int i = 0; i < image_size[0]*image_size[1]; i++){

            image[i] *= (region[i] == 0);
        }
        write_bmp(image, image_size[0], image_size[1]);
    }
}


int main(int argc, char** argv){
    
    init_mpi(argc, argv);
    
    load_and_allocate_images(argc, argv);
    
    create_types();
    
    distribute_image();

    grow_region();
    
    gather_region();
    
    MPI_Finalize();
    
    write_image();
    
    exit(0);
}
