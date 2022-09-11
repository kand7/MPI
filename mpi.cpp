//VIKTOR ROMANIOYK//
//713242017024//
//TMHMA E8//
//COMPILE ME -lm (gia na doylepsei to POW)//
//p.x. mpicc mpi_ask_2.cpp -o mpi_ask_2 -lm//
 



#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
int main(int argc, char** argv) 
{
	MPI_Init(&argc,&argv);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int p; 
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	int num_elements;
	int num_elements_per_proc;
	int array[100];
	int total_max;                     //DHLWSH VASIKWN METABLITWN//
					   //POY THA XREISIMOPOIHTHOYN//
											
	int total_min=0;
	int k;
	int max_counter=0;
	int min_counter=0;
	if (my_rank== 0)  
	{	printf("Dose plithosarithmvn:\n"); 
		scanf("%d", &num_elements);
		printf("Dose ta stoixeia \n");
		int i;				//DWSE STOIXEIA//
		for (i=0;i<num_elements;i++)
		{
			scanf("%d", &array[i]);
		}
	}
	//STEILE STOYS ALLOYS TA VASIKA STOIXEIA//
	MPI_Bcast(&num_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&array, 1, MPI_INT, 0, MPI_COMM_WORLD);
	num_elements_per_proc= num_elements/p;		
	//*SUB_ARRAY -> BUFFER//
	int *sub_array= (int *)malloc(sizeof(int) * num_elements_per_proc);
	int *prefix_array= (int *)malloc(sizeof(int) * num_elements_per_proc);
	//KANE SCATTER//
	MPI_Scatter(array, num_elements_per_proc, MPI_INT, sub_array, num_elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	int sum = 0;
	int sub_avg;
	int l;
	for (l = 0; l < num_elements_per_proc; l++)
	{
		sum += sub_array[l];
	}
	sub_avg= sum/ (int)num_elements_per_proc;
	int *sub_avgs;
	if (my_rank== 0) 
	{
		sub_avgs= (int *)malloc(sizeof(int) * p);
	}
	//KANE TO GATHER//
	MPI_Gather(&sub_avg, 1, MPI_INT, sub_avgs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//KANE SCATTER//
	MPI_Scatter(array, num_elements_per_proc, MPI_INT, prefix_array, num_elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	//PREFIXSUM//
	//ALLA EMFANIZEI APLA TA STOIXEIA, ERROR//
	int prefixsum[num_elements_per_proc];
	prefixsum[0]=prefix_array[0];
	int q;
	for (q = 1; q < num_elements_per_proc; q++)
	{
		prefixsum[q]=prefixsum[q-1]+prefix_array[q];
	}
	int *prefixarray;
	if (my_rank== 0) 
	{
		prefixarray= (int *)malloc(sizeof(int) * p);

	}
	//GATHER//
	MPI_Gather(&prefixsum, 1, MPI_INT,prefixarray, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//STEILE TA DEDOMENA SE OLOYS//
	MPI_Bcast(&num_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&array, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int local_max=0;
	for (k=0; k<num_elements; k++)
	{
		if (array[k]>local_max) local_max=array[k];
	}

	//SYLOGH ME KRITHRIO TO MAX//
	MPI_Reduce(&local_max, &total_max, 1, MPI_INT,MPI_MAX,0, MPI_COMM_WORLD);
	//STEILE TA DEDOMENA SE OLOYS//
	MPI_Bcast(&num_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&array, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int local_min=999999999;
	for (k=0; k<num_elements; k++)
	{
		if (array[k]<=local_min) local_min=array[k];
	}
	//SYLOGH ME KRIHRIO TO MIN(MAX)//
	//GIA KAPOIO LOGO AN VALW MPI_MIN DEN MOY DOYLEVEI SWSTA????//
	MPI_Reduce(&local_min, &total_min, 1, MPI_INT,MPI_MAX,0, MPI_COMM_WORLD);
	//KWDIKAS ARXIGOY//
	if (my_rank== 0)
	{
		float avg;
		int tot_sum= 0;
		int j;
		int var_pow=0;   //DHLWSH METAVLITWN//
		float var_total=0;
		float diany=0;
		float diany_max=0;
		int diany_num=0;
		for (j = 0; j < p; j++) 
		{
			tot_sum+= sub_avgs[j]; //YPOLOGISMOS TOY SUM OF ELEMENTS//
		}
		avg= (float)tot_sum/ (float) p; //TO M//
		printf("Avgof all elements is %f\n", avg);
		printf("SUM of Elements is %d\n", tot_sum);
		int i;
		for (i=0;i<num_elements;i++)
		{
		if (array[i]<avg) min_counter+=1; //POSA MIKROTERA APO TO M//
		if (array[i]>avg) max_counter+=1; //POSA MAGALYTERA APO TO M//
		}
		for (i=0;i<num_elements;i++)
		{
		var_pow+= pow( (array[i] - avg), 2); //YPOLOGISMOS TOY VAR//
		}
		var_total=float(var_pow)/float(num_elements);
		//YPOLOGISMOS TOY DIANISMATOS//
		for (i=0;i<num_elements;i++)
		{
		diany=( float((array[i] - total_min))/float((total_max-total_min))*100);
		printf("To DIANISMA %d EINAI %f\n",i,diany);
		//POIO EXEI TO MEGALYTERO DIANYSMA//
		if (diany>diany_max)
		{
			diany_max=diany;
			diany_num=array[i];
			
		}
		}
		//EKTYPWSEIS//
		printf(" VAR IS : %f \n",var_total);
		
		printf(" %d Stoixeia einai mikrotera apo to m \n",min_counter);
		printf(" %d Stoixeia einai megalytera apo to m \n",max_counter);	
		
		printf("MAX IS %d\n", total_max);
		printf("Min IS %d\n", total_min);
		printf("MEGALYTERO DIANYSMA EINAI %f TOY STOIXEIOY %d\n",diany_max,diany_num);
		//PREFIXSUM TO OPOIO DEN DOYLEVEI//
		for (i=0;i<num_elements;i++)
		{
		printf("PREFIX ARRAY %d = %d \n",i,prefixarray[i]);
		}
	}
		
	MPI_Finalize();
	return 0;
	//TELOS//
}
