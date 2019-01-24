#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <mpi.h>

//Compile with mpicc run_clones_MPI.c -o run_clones_MPI
//Run with mpirun -np # run_clones_MPI


int main(){
  MPI_Init(NULL, NULL);

  float HII_EFF_FACTOR[] = {20, 27.79, 38.61, 53.65, 74.55, 103.59, 143.94, 200};
  float R_BUBBLE_MAX[] = {5, 6.60, 8.72, 11.51, 15.20, 20.07, 26.51, 35};
  float ION_Tvir_MIN[] ={8e3, 1.15e4, 1.65e4, 2.36e4, 3.39e4, 4.86e4, 6.97e4, 1e5};
  //Uncomment for the 8x8x8 sampling

  /*//Uncomment for the 10x10x10 sampling
  float HII_EFF_FACTOR[] = {10, 13.95, 19.45, 25.14, 37.14, 52.82, 73.68, 102.78, 143.37, 200};
  float R_BUBBLE_MAX[] = {10, 12.51, 15.65, 19.57, 24.48, 30.63, 38.32, 47.93, 59.96, 75};
  float ION_Tvir_MIN[] ={8e3, 1.06e4, 1,40e4, 1.86e4, 2.46e4, 3.25e4, 4.31e4, 5.70e4, 7.55e4, 1e5};
  */
  /*//Uncomment for the 20x6x20 sampling
  float HII_EFF_FACTOR[] = {20, 22.58, 25.49, 28.77, 32.48, 36.66, 41.33, 46.71, 52.73, 59.53, 67.20, 75.85, 85.63, 96.66, 109.11, 123.17, 139.04, 156.95, 177.17, 200};
  float R_BUBBLE_MAX[] = {5, 7.38, 10.89, 16.07, 23.72, 35};
  float ION_Tvir_MIN[] ={8e3, 9.14e3, 1.04e4, 1.19e4, 1.36e4, 1.56e4, 1.78e4, 2.03e4, 2.32e4, 2.65e4, 3.02e4, 3.45e4, 3.94e4, 4.50e4, 5.14e4, 5.88e4, 6.71e4, 7.67e4, 8.76e4, 1e5};
  */
  //based on Mesinger & Grieg 2017
  char cmnd[1000];
  int i,j,l,i_max,j_max,l_max,threads,id,clone,total_clones;
  int nb_tasks, my_rank, err_mpi, nb_lines, line;
  float v1,v2,v3,w1,w2,w3;
  FILE *file;

  //USE THIS TO DECIDE HOW THE CODE WORKS!!
  int NORMAL = 1; //Uses the parameter space defined above
  int RANDOM_CHOICE = 0; //Randomly choose parameter values
  int READ_IN = 0; //Reads in parameter values

  err_mpi = MPI_Comm_size(MPI_COMM_WORLD, &nb_tasks);
  err_mpi = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  i_max = sizeof(HII_EFF_FACTOR)/sizeof(HII_EFF_FACTOR[0]);
  j_max = sizeof(R_BUBBLE_MAX)/sizeof(R_BUBBLE_MAX[0]);
  l_max = sizeof(ION_Tvir_MIN)/sizeof(ION_Tvir_MIN[0]);

  total_clones = i_max*j_max*l_max;
  if (RANDOM_CHOICE) total_clones = nb_tasks;

  clone = my_rank;
  i = clone % i_max;
  j = clone / i_max % j_max;
  l = clone / (i_max*j_max);

  if (RANDOM_CHOICE){
  	srand(my_rank*1000);
	HII_EFF_FACTOR[i] = ((double) rand() / (double) RAND_MAX)*(HII_EFF_FACTOR[i_max-1] - HII_EFF_FACTOR[0]) + HII_EFF_FACTOR[0];
	R_BUBBLE_MAX[j] = ((double) rand() / (double) RAND_MAX)*(R_BUBBLE_MAX[j_max-1] - R_BUBBLE_MAX[0]) + R_BUBBLE_MAX[0];
	ION_Tvir_MIN[l] = ((double) rand() / (double) RAND_MAX)*(ION_Tvir_MIN[l_max-1] - ION_Tvir_MIN[0]) + ION_Tvir_MIN[0];
  }

  //If the file of resampled points is to be read:
  nb_lines = 1;
  if (READ_IN){
        nb_lines--;
  	file = fopen("resampled.dat","r");
	while (fscanf(file, "%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n",&i,&j,&l,&v1,&v2,&v3,&w1,&w2,&w3) == 9){
		if (w1 != 0 || w2 != 0 || w3 != 0){nb_lines++;}
  	}
  	fclose (file);
  }
  float resampled[nb_lines][6];
  if (READ_IN){
  	line = 0;
  	file = fopen("resampled.dat","r");
	while (fscanf(file, "%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n",&i,&j,&l,&v1,&v2,&v3,&w1,&w2,&w3) == 9){
		if (w1 != 0 || w2 != 0 || w3 != 0){
			resampled[line][0] = (float) i;
			resampled[line][1] = (float) j;
  			resampled[line][2] = (float) l;
			resampled[line][3] = w1;
			resampled[line][4] = w2;
  			resampled[line][5] = w3;
			line++;
		}
  	}
  	fclose (file);
  } else {float resampled[1][6];}
 
  //Make sure the clone directory is empty
  //system("rm -rf 21cmFAST-clones/*");

  if (NORMAL) fprintf(stderr, "Creating clone %i/%i with parameters f = %.0f, R = %.0f, T = %.1e\n",clone, total_clones, HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  if (RANDOM_CHOICE) fprintf(stderr, "Creating clone %i/%i with parameters f = %.0f, R = %.0f, T = %.1e\n",clone, nb_tasks, HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  if (READ_IN) fprintf(stderr, "Creating clone %i/%i with parameters f = %.0f, R = %.0f, T = %.1e\n",clone, nb_lines, resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  //Copy 21cmFAST-master into the clones directory and rename the clone with the relative parameters
  if (!READ_IN) sprintf(cmnd,"cp -R 21cmFAST-master 21cmFAST-clones/f%.0f_R%.0f_T%.1e", HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  if (READ_IN) sprintf(cmnd,"cp -R 21cmFAST-master 21cmFAST-clones/f%.0f_R%.0f_T%.1e", resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  system(cmnd);
  //Update the parameter files (The values in the master should be 15, 30, 1e4 respectively for this to work)
  if (!READ_IN){
  	sprintf(cmnd, "sed -i 's/HII_EFF_FACTOR (float) (30)/HII_EFF_FACTOR (float) (%.2f)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", HII_EFF_FACTOR[i], HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  	system(cmnd);
  	sprintf(cmnd, "sed -i 's/R_BUBBLE_MAX (float) (50)/R_BUBBLE_MAX (float) (%.2f)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", R_BUBBLE_MAX[j], HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  	system(cmnd);
  	sprintf(cmnd, "sed -i 's/ION_Tvir_MIN (double) (3e4)/ION_Tvir_MIN (double) (%.2e)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", ION_Tvir_MIN[l], HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  	system(cmnd);
  }else{
  	sprintf(cmnd, "sed -i 's/HII_EFF_FACTOR (float) (30)/HII_EFF_FACTOR (float) (%.2f)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", resampled[clone][3], resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  	system(cmnd);
  	sprintf(cmnd, "sed -i 's/R_BUBBLE_MAX (float) (50)/R_BUBBLE_MAX (float) (%.2f)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", resampled[clone][4], resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  	system(cmnd);
  	sprintf(cmnd, "sed -i 's/ION_Tvir_MIN (double) (3e4)/ION_Tvir_MIN (double) (%.2e)/g' 21cmFAST-clones/f%.0f_R%.0f_T%.1e/Parameter_files/ANAL_PARAMS.H", resampled[clone][5], resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  	system(cmnd);
  }
  //Go into the new clone's Program directory and run 'make'
  if (!READ_IN) sprintf(cmnd,"21cmFAST-clones/f%.0f_R%.0f_T%.1e/Programs", HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  if (READ_IN) sprintf(cmnd,"21cmFAST-clones/f%.0f_R%.0f_T%.1e/Programs", resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  chdir(cmnd);
  system("make");
  //system("rm *.c");
  //system("rm PROGRAM_LIST");
  //system("rm Makefile");
  //chdir("../../..");

  if (NORMAL) fprintf(stderr, "Treating clone = %i/%i with thread %i\n", clone,total_clones,clone);
  if (RANDOM_CHOICE) fprintf(stderr, "Treating clone = %i/%i with thread %i\n", clone,nb_tasks,clone);
  if (READ_IN) fprintf(stderr, "Treating clone = %i/%i with thread %i\n", clone,nb_lines,clone);

  //Go to the clone directory
  //if (!READ_IN) sprintf(cmnd,"21cmFAST-clones/f%.0f_R%.0f_T%.1e/Programs", HII_EFF_FACTOR[i], R_BUBBLE_MAX[j], ION_Tvir_MIN[l]);
  //if (READ_IN) sprintf(cmnd,"21cmFAST-clones/f%.0f_R%.0f_T%.1e/Programs", resampled[clone][3], resampled[clone][4], resampled[clone][5]);
  //chdir(cmnd);
  //Run 21cmFAST
  system("./drive_zscroll_noTs");
  chdir("..");
  system("mv Output_files/Deldel_T_power_spec/* .");
  //system("rm -r */");
  chdir("../..");

  //fprintf(stderr, "Clone %i/%i complete\n *****************************\n", clone, total_clones);
  MPI_Finalize();

  return 0;
}
