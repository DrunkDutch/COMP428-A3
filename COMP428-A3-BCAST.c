/***
COMP428 Assignment 3
Devin Mens 26290515
**/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#define MASTER 0        /* task ID of master task */

// http://stackoverflow.com/questions/3437404/min-and-max-in-c
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

//Flatten 2D matrix coords into 1D array offset
int getIndex(int x, int y, int n)
{
    return x + (y * n);
}

//Splits a 1D array mapping into the proper 2D array/matrix
void getCoords(int i, int n, int *point)
{
    point[0] = i % n;
    point[1] = i / n;
}

//Apparently something is giving me negative numbers, possible int overflow
int maxIntAdd(int a, int b)
{
    if(a == INT_MAX || b == INT_MAX)
        return INT_MAX;

    return a + b;
}

// Helper to clean up code for array instantiation (I hate malloc)
int* gen_Array(int size){
    return malloc(size * sizeof(int));
}

int main (int argc, char *argv[])
{
int	taskid,	        /* task ID - also used as seed number */
	numtasks,       /* number of tasks */
	rc,             /* return code */
	i;

    struct timeval start, end;

    // Record the start time
    gettimeofday(&start, NULL);

    MPI_Status status;

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);


    int pn = sqrt(numtasks);
    if(!((pn * pn) == numtasks))
    {
        if(taskid == MASTER)
            printf ("The number of processes must be a square number\n");

        MPI_Abort(MPI_COMM_WORLD, 1);
        return rc;
    }

    int* inputValue;
    int n;
    int inputSize = 0;

    // Only the master will deal with reading the file
    if(taskid == MASTER)
    {
        // Open the input file
        FILE *myFile;
        myFile = fopen("input.txt", "r");

        int input;
        while(!feof(myFile))
        {
            int c = fscanf(myFile, "%d\t", &input);

            if(c != 1)
            {
                char word[16];
                fscanf(myFile, "%s\t", &word);
            }

            inputSize++;
        }

        n = sqrt(inputSize);
        inputValue = gen_Array(inputSize);

        rewind(myFile);

        int count = 0;
        while(!feof(myFile))
        {
            int c = fscanf(myFile, "%d", &inputValue[count]);

            if(c != 1)
            {
                inputValue[count] = INT_MAX;
                char word[16];
                fscanf(myFile, "%s", &word);
            }

            count++;
        }

        fclose(myFile);
    }

    // Broadcast the number of inputs
    MPI_Bcast(&inputSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    if(pn > n)
    {

        printf ("Currently only works with 1 node / 1 processor mapping\n");

        MPI_Abort(MPI_COMM_WORLD, 1);
        return rc;
    }

    // Non-master processes will need to create the sequence buffer for themselves
    if(taskid != MASTER)
    {
        inputValue = gen_Array(inputSize);
    }

    // Broadcast the sequence to all processors
    MPI_Bcast(inputValue, inputSize, MPI_INT, MASTER, MPI_COMM_WORLD);


    int subsection_size[pn];    // The number of rows and columns a section along the axis responsible for
    int offset[pn];     // The x and y coordinate where a row or column section starts

    offset[0] = 0;

    // Calculate the sections of the table each process will be responsible for
    int baseCount = n / pn;
    int remainder = n % pn;

    for(i = 0; i < pn; i++)
    {
        subsection_size[i] = baseCount;

        if(i < remainder)
        {
            subsection_size[i]++;
        }

        if(i != 0)
        {
            offset[i] = offset[i - 1] + subsection_size[i - 1];
        }
    }

    // Find the corresponding k for a processor's dimension
    int kcounts[n];
    int ncount = 0;
    for(i = 0; i < n; i++)
    {
        if(offset[ncount] + subsection_size[ncount] == i)
            ncount++;

        kcounts[i] = ncount;
    }

    int point[2];
    getCoords(taskid, pn, point);

    // Place each process into two comm worlds, for row and column
    // Because MPI_Cart_sub doesnt actually work
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, point[1], point[0], &row_comm);

    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, point[0], point[1], &col_comm);

    int k, x, y;
    for(k = 0; k < n; k++)
    {
        // Find the processor dimension responsible for this k
        int dim_needed = kcounts[k];

        // Allocate the buffers
        int* rowbuffer = gen_Array(subsection_size[point[0]]);
        int* colbuffer = gen_Array(subsection_size[point[1]]);

        //Processor holds the k row needed
        //Generate buffer of row values to be sent along column
        if(point[1] == dim_needed)
        {
            for(i = 0; i < subsection_size[point[0]]; i++)
            {
                int index = getIndex(offset[point[0]] + i, k, n);
                rowbuffer[i] = inputValue[index];
            }
        }

        // Perform the column wise broadcast of the row
        MPI_Bcast(rowbuffer, subsection_size[point[0]], MPI_INT, dim_needed, col_comm);

        // Same as with row but with column
        if(point[0] == dim_needed)
        {
            for(i = 0; i < subsection_size[point[1]]; i++)
            {
                int index = getIndex(k, offset[point[1]] + i, n);
                colbuffer[i] = inputValue[index];
            }
        }

        // Perform the row wise broadcast of the column
        MPI_Bcast(colbuffer, subsection_size[point[1]], MPI_INT, dim_needed, row_comm);

        // Locally perform the Floyd's all pair calculation
        for(y = 0; y < subsection_size[point[1]]; y++)
        {
            for(x = 0; x < subsection_size[point[0]]; x++)
            {
                // Calculate actual matrix coords from 1D array
                int m_x = x + offset[point[0]];
                int m_y = y + offset[point[1]];
                if(m_x == k || m_y == k || m_x == m_y)
                    continue;

                inputValue[getIndex(m_x, m_y, n)] =
                    min(inputValue[getIndex(m_x, m_y, n)],
                    maxIntAdd(rowbuffer[x], colbuffer[y]));
            }
        }

    }

    // The master gathers all data to generate new matrix
    if(taskid == MASTER)
    {
        for(i = 1; i < numtasks; i++)
        {
            // Find the virtual locations of processor
            int vir_coords[2];
            getCoords(i, pn, vir_coords);

            int receiveBufferSize = subsection_size[vir_coords[0]] * subsection_size[vir_coords[1]];
            int* receiveBuffer = gen_Array(receiveBufferSize);

            MPI_Recv(receiveBuffer, receiveBufferSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(y = 0; y < subsection_size[vir_coords[1]]; y++)
            {
                for(x = 0; x < subsection_size[vir_coords[0]]; x++)
                {
                    inputValue[getIndex(offset[vir_coords[0]] + x, offset[vir_coords[1]] + y, n)] = receiveBuffer[x + (y * subsection_size[vir_coords[0]])];
                }
            }

        }

        FILE *myFile;
        myFile = fopen("output.txt", "w");

        // Write the result to output.txt
        int x,y;
        for(y = 0; y < n; y++)
        {
            for(x = 0; x < n; x++)
            {
                fprintf(myFile, "%d\t", inputValue[getIndex(x, y, n)]);
            }

            fprintf(myFile, "\n");
        }

        fclose(myFile);
    }
    else
    {
        // Slaves determine responsibility and send to master as needed
        // Moved this loop to after all k matrices are calculated. Getting weird errors
        int sendBufferSize = subsection_size[point[0]] * subsection_size[point[1]];
        int* sendBuffer = gen_Array(sendBufferSize);

        for(y = 0; y < subsection_size[point[1]]; y++)
        {
            for(x = 0; x < subsection_size[point[0]]; x++)
            {
                sendBuffer[x + (y * subsection_size[point[0]])] = inputValue[getIndex(offset[point[0]] + x,
                                                                              offset[point[1]] + y, n)];
            }
        }

        MPI_Send(sendBuffer, sendBufferSize, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

    }

    MPI_Finalize();

    if(taskid == MASTER)
    {
        // Temp code for console out

        printf("The final buffer is:\n");
        for(k = 0; k < n; k++)
        {
            for(i = 0; i < n; i++)
            {
                printf("%11d", inputValue[getIndex(i, k, n)]);
            }
            printf("\n");
        }


        gettimeofday(&end, NULL);

        printf("took %lu seconds and \n", end.tv_sec - start.tv_sec);
        printf("took %lu microseconds\n", end.tv_usec - start.tv_usec);

    }


    return 0;
}