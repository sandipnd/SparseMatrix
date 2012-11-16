#include<iostream>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include <time.h>
#include <stddef.h>
#include <sys/time.h>
#include <unistd.h>  
#include <fstream>

#define MAX(a,b) ( a>b?a:b)
#define MIN(a,b) (a<b?a:b)
#define PI 3.14159265
#define VAL .01
#define WHERE fprintf(stdout,"In line no = %d\n",__LINE__);
#define REPITITIONS 100


double *xVector , *yVector;
using namespace std;

                                                                
 
 double elapsed_time()
 {
    struct timeval t;
     struct timezone whocares;
     double total;
     double sec, msec;     /* seconds */
     double usec;          /* microseconds */
 
//     // gettimeofday(&t, NULL);
     gettimeofday(&t, &whocares);
// 
     msec = (double) (t.tv_sec);
     usec = 1.0e-6*(double) (t.tv_usec);
     total = msec + usec;
     if (total < 0) 
         return(-17.0);
     else
         return(total);
 
 }


void COOMultiply(int nnz ,int n, int *row, int *column , double *value ) {
  
  for(int i = 1; i <= n ; i++ ) 
     yVector[i] =0;

  for(int k = 1; k <= nnz ; k ++ ) {
     // cout << value[k] << ":" << row[k]<< ":" << column[k] << ":" << k << endl;
      yVector[row[k]] = yVector[row[k]] + value[k] * xVector[column[k]]; 
  }
 /*
  cout << "\n\n\n";
  for(int k = 1; k <= n; k ++ ) 
   cout << yVector[k] << endl;
 cout << "\n";
 */
}

void CSRMultiply(int nnz ,
			int n, 
			  int *row, 
			    int *column , 
				double *value ) {
 // cout << "********************************\n";  
  for(int i = 1; i <= n ; i++ ) 
     yVector[i] = 0;
    
  for(int i = 1; i <= n ; i++ ) {
     for ( int k = row[i]; k <= row[i+1] - 1; k++ ) {
	//cout << "value :"  << value[k] << "\t";
       yVector[i] = yVector[i] + value[k] * xVector[column[k]];
    } 
  }
  /*
  cout << "------------------------------------------\n";
  for(int k = 1; k <= n; k ++ ) 
   cout << yVector[k] << endl;
   */
}
void CSCMultiply(int nnz , int n ,int *row, int *column , double *value ) {
  
 for(int i = 1; i <= n ; i++ ) 
     yVector[i] = 0;
  
 for(int i = 1; i <= n ; i++) {
    for(int j = column[i]; j <= column[i+1] - 1; j++) {
       yVector[row[j]] = yVector[row[j]] + value[j] * xVector[i];
    }
  }
  /*cout << "---------------------CSC---------------------\n";
  for(int k = 0; k <= n; k ++ ) 
   cout << yVector[k] << endl;
  */
}

void coocsr(int **newrow , int **newcolumn , double **newvalue ,int *oldrow, int *oldcolumn , double *oldvalue , int nnz , int n) {
  
  
  //for(int k = 1; k <= nnz ; k++ )
   //cout << oldrow[k] << "\t";

  //cout << "\n---------------------\n";
  //WHERE
  int *crow = *newrow;
  for ( int  i =1; i <= n; i++ ) {
     crow[i] = 0;
  }
  for(int j = 1; j <= nnz ; j++ ) { 
    crow[oldrow[j]] =  crow[oldrow[j]] + 1;
  }
  //WHERE
  int tmpVal = 1;
  for(int k = 1 ; k <= n; k++) {
   int k0 =  crow[k];
   crow[k] = tmpVal;
   tmpVal += k0;
   //crow[k+1] = crow[k];
  }
  //WHERE
   // cout << "nnzs is = " << nnz << "\t";
  ////WHERE
  int tmp, iad;
 for(int k = 1 ; k < nnz; k++) {
     tmp = oldrow[k];
     int iad = crow[tmp];
    (*newvalue)[iad]= oldvalue[k];
    (*newcolumn)[iad] = oldcolumn[k];
    crow[tmp] = iad + 1;
  }
  //WHERE
  for(int j = n; j >= 1; j-- )
    crow[j+1] = crow[j];
  
  crow[1]= 1;
  /*
  cout << "\n";
  for(int i = 1 ; i <= nnz ; i++)
   cout << (*newvalue)[i] << "\t";

  cout << "---------------\n";
  */
}
void csrcsc(int *iao, int *jao, double *ao,int *ia, int *ja, double *a,int nnz , int n) {

   
  for(int i = 1; i <= n ; i++)
    iao[i] = 0;

  for (int i = 1; i <= n ; i++ ) {
   for(int k = ia[i]; k <=  ia[i+1] - 1 ; k++) { 
      int j   =  ja[k] + 1;
      iao[j] =  iao[j] + 1;
      }
  }
  iao[1]  =  1;
  for (int i= 1; i <= n ; i++ )
    iao[i+1]  =  iao[i] + iao[i+1];
    
 for(int i = 1; i <= n ; i ++ ) {
   for(int k = ia[i];k <= ia[i+1] - 1 ; k++) {
     int j        =  ja[k]; 
     int next     =  iao[j]; 
     //cout <<  next <<  " : " << k << endl;
     
     ao[next]     =  a[k];   
     //cout << "val : " << ao[next] << endl; 
     jao[next]    =  i;
     iao[j]       =  next + 1;
     }
  }
  
 for( int i = n; i >= 1 ; i-- ) {  
   iao[i+1] = iao[i];
 }
 ////WHERE
 iao[1]  =  1; 
 

/*
 for(int i = 1; i <=nnz ; i++ )
   cout << ao[i] << "\t";
*/
}

int main(int argc , char *argv[]) {

 int b = 5, halfb, nzcount = 0;
 int i , j , index = 1 , n = 100 , nhigh = 10000, nhzcount = 0;

 int *csrcscrow = NULL, *csrcsccolumn = NULL;
 double  *csrcscvalue = NULL; 
 double starttime, endtime, looptime;
 double timeDiff1,timeDiff2;
// char *filename = "data_overheads";
 double gflops_per_sec , cost_in_flops , timer_overhead;
 FILE *fp;
 fp = fopen("logtime.txt","a");
 ifstream myfile ("numberRange");
  string init , incr , last ;
  if (myfile.is_open()) {
    getline (myfile,init);
    getline (myfile,incr);
    getline (myfile,last);
    myfile.close();
 }
 
 nhigh = atoi(last.c_str());
 yVector = new double[nhigh+1];
 xVector = new double[nhigh+1];
 
 
  halfb = floor(b/2);
  for(i = 1; i <= nhigh ; i++) {
    for(j= MAX(1,(i-halfb)); j <= MIN(nhigh,(i+halfb)); j++ ) {
     nhzcount += 1; 
  }
 }
  csrcscrow = new int[nhigh+1];
  csrcsccolumn = new int[nhzcount+1];
  csrcscvalue = new double[nhzcount+1];  
  int *iao = new int[nhzcount+1];
  int *jao = new int[nhigh+1];
  double *ao = new double[nhzcount+1];

  //cout << nzcount;
  int *row, *column ;
  double  *value;
  //cout << "nhz " << nhzcount << "\t";
  row = new  int[nhzcount+1];
  column = new int[nhzcount+1];
  value = new double[nhzcount+1];

  
  //cout << init << " " << incr << " " << last;
  
  
 //start of loop for n 
 for(int k = 1; k <= nhigh ; k ++ ) 
   xVector[k] = VAL; 

for ( int loop = atoi(init.c_str()); loop < nhigh ; loop += atoi(incr.c_str()))  {

 n = loop;
 nzcount = 0;
 index = 1;
  for ( int rst = 1; rst <= nhzcount ; rst ++ ) { 
    row[rst] = 0;
    column[rst] = 0;
    value[rst] = 0.0;
    csrcscvalue[rst] = 0.0;
    csrcsccolumn[rst] = 0;
    iao[rst] = 0;
     ao[rst] =0;
   }
  for ( int rst = 1; rst <= nhigh ; rst ++ ) {
      csrcscrow[rst] = 0;
      jao[rst] =0; 
  }
 for(i = 1; i <= n ; i++) {
    for(j= MAX(1,(i-halfb)); j <= MIN(n,(i+halfb)); j++ ) {
     nzcount += 1; 
  }
 }

 //cout << " COO Matrix formation started\n";
 for(i = 1; i <= n ; i++) {
   for(j= MAX(1,(i-halfb)); j <= MIN(n,(i+halfb)); j++ ) {
     //cout << "i " << i << " : " << j << endl;
      row[index] = i;
      column[index] = j;
      value[index] = ( PI + i ) / (double)j;
      index++;   
   }
  }
 //cout << " COO Matrix formation End\n";
   //cout << " COO -Vector Multiplication starting for n = : " << n << endl ;
   starttime = elapsed_time();   
   for(int count = 0 ; count < REPITITIONS ; count++)
   COOMultiply( nzcount,n,row,column,value);
   endtime =  elapsed_time();
   timeDiff1 =  endtime - starttime;
   gflops_per_sec=  ((2.0e-9)  * n * REPITITIONS) / timeDiff1; 
   double avgtime = timeDiff1 / REPITITIONS;
   fprintf(fp,"COO : n = %d : time = %.17e : gflops = %.17e \n",n,avgtime,gflops_per_sec);
   //cout << " COO -Vector Multiplication End\n";

   //cout << " COO to CSR conversion starting\n";
   ////WHERE
   coocsr(&csrcscrow,&csrcsccolumn,&csrcscvalue,row,column,value,nzcount,n);
   //cout << " COO to CSR conversion Ended\n";
   

   //cout << " CSR -Vector Multiplication starting for n= : " << n << endl;
   starttime = elapsed_time(); 
   for(int count = 0 ; count < REPITITIONS ; count++)
   CSRMultiply(nzcount,n,csrcscrow,csrcsccolumn,csrcscvalue);   
   endtime =  elapsed_time();
   timeDiff1 =  endtime - starttime;
   gflops_per_sec=  ((2.0e-9)  * n * REPITITIONS) / timeDiff1; 
   avgtime = timeDiff1 / REPITITIONS;
   fprintf(fp,"CSR : n = %d : time = %.17e : gflops = %.17e \n",n,avgtime,gflops_per_sec);
   //cout << " CSR -Vector Multiplication End\n";
   //cout << "nzcount " << nzcount << "\t";
   
   //double ao[nzcount+1];
   //   cout << " CSR to CSC conversion starting\n";
   csrcsc(jao,iao,ao,csrcscrow,csrcsccolumn,csrcscvalue,nzcount,n);
   //cout << " CSR to CSC conversion  End\n";
   
   //cout << " CSC -Vector Multiplication starting for n= : " << n << endl;
   starttime = elapsed_time(); 
   for(int count = 0 ; count < REPITITIONS ; count++)
   CSCMultiply(nzcount,n,iao,jao,ao); 
   //WHERE
   endtime =  elapsed_time();
   timeDiff1 =  endtime - starttime;
   gflops_per_sec=  ((2.0e-9)  * n * REPITITIONS) / timeDiff1; 
   avgtime = timeDiff1 / REPITITIONS;
   fprintf(fp,"CSC : n = %d : time = %.17e : gflops = %.17e \n",n,avgtime,gflops_per_sec);
   //cout << " CSC -Vector Multiplication End\n";
   fprintf(fp,"-----------------------------------------------\n");   
   //cout << "----------End of loop-------------- for n = " << loop << endl;
   
   }
   //WHERE  
   delete[] iao;
   delete[] jao;
   delete[] ao;
   delete[] yVector;
   delete[] xVector;
   delete[] row;
   delete[] column;
   delete[] value;
   
   delete[] csrcscrow;
   delete[] csrcsccolumn;
   delete[] csrcscvalue;
   fclose(fp);  
}


