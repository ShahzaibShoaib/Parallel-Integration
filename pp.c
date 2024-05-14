#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "mpi.h"

#define NARGS 4
#define Z0 0.000000001

/* Function prototypes */
double select_integral_method(int rule_flag, int fun_flag, double local_a, double local_b, int local_n, double h);
void Get_data(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p);
double Rectangle(double local_a, int local_n, double h, double (*f)(double x));
double Trapezoidal(double local_a, double local_b, int local_n, double h, double (*f)(double x));
double Simpsons(double local_a, int local_n, double h, double (*f)(double x));
double Rectangle2(double local_a, int local_n, double h, double (*f)(double x));
double Trapezoidal2(double local_a, double local_b, int local_n, double h, double (*f)(double x));
double pi(double x);
double ex2(double x);
double sinx(double x);
double x2(double x);

int main(int argc, char** argv) {
    int my_rank, p, n, local_n, i, source;
    int dest = 0, tag = 0, fun_flag = 0, rule_flag = 0;
    double a, b, h, local_a, local_b, integral, total;
    double overhead, start, finish, cpu_precision;
    int ch;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != NARGS) {
        rule_flag = 2;
        fun_flag = 2;
        n = 100000;
        printf("Note: Incorrect parameter(s). Running with default rule.\n");
        printf("Usage: %s [-rule -function N]\n", argv[0]);
        printf("\nFunctions: Rules:\n");
        printf("********* ******\n");
        printf(" -p = pi -r = Rectangle\n");
        printf(" -e = ex2 -s = Trapezoidal\n");
        printf(" -i = sinx -t = Simpsons\n");
        printf(" -x = x2 -o = Optimized rectangle\n");
        printf(" -m = Modified trapezoidal\n");
        printf("N = number of segments (rectangles, trapezoids or arcs) \n");
    } else {
        n = atol(argv[NARGS-1]);
        while ((ch = getopt(argc, argv, "piestormx")) != -1) {
            switch (ch) {
                case 'r': rule_flag = 1; break;
                case 't': rule_flag = 2; break;
                case 's': rule_flag = 3; break;
                case 'o': rule_flag = 4; break;
                case 'm': rule_flag = 5; break;
                case 'p': fun_flag = 1; break;
                case 'e': fun_flag = 2; break;
                case 'i': fun_flag = 3; break;
                case 'x': fun_flag = 4; break;
                default: printf("Usage: %s [-rule -function]\n", argv[0]); return(2);
            }
        }
    }

    a = Z0;
    b = 1;
    cpu_precision = MPI_Wtick();

    overhead = 0.0;
    for (i = 0; i < 100; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        finish = MPI_Wtime();
        overhead += (finish - start);
    }
    overhead /= 100.0;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    Get_data(&a, &b, &n, my_rank, p);
    h = (b-a)/n;
    local_n = n/p;
    local_a = a + my_rank * local_n * h;
    local_b = local_a + local_n * h;

    // Integration decision tree
    integral = select_integral_method(rule_flag, fun_flag, local_a, local_b, local_n, h);

    if (my_rank == 0) {
        total = integral;
        for (source = 1; source < p; source++) {
            MPI_Recv(&integral, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
            total += integral;
        }
    } else {
        MPI_Send(&integral, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();

    if (my_rank == 0) {
        printf("\nIntegration result (area) = %.25f\n", total);
        printf("Start time: %.10f, Finish time: %.10f, difference: %.10f, overhead: %.10f, precision is %.10f\n",
               start, finish, (finish - start), overhead, cpu_precision);
        printf("==========================================================\n");
    }

    MPI_Finalize();
    return 0;
}

// Function implementations would follow here...

/* Function that distributes data to all processes */
void Get_data(double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p) {
    int source = 0; // Source process for MPI
    int dest; // Destination process for MPI
    int tag;
    MPI_Status status; // Declare MPI_Status variable

    if (my_rank == 0) {
        // Process 0 sends data to all other processes
        for (dest = 1; dest < p; dest++) {
            tag = 0;
            MPI_Send(a_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            tag = 1;
            MPI_Send(b_ptr, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
            tag = 2;
            MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
        }
    } else {
        // All other processes receive data from process 0
        tag = 0;
        MPI_Recv(a_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        tag = 1;
        MPI_Recv(b_ptr, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        tag = 2;
        MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
    }
}

/* Rectangle rule integration */
double Rectangle(double local_a, int local_n, double h, double (*f)(double x)) {
    double integral = 0.0;
    double x;
    for (int i = 0; i < local_n; i++) {
        x = local_a + i * h;
        integral += f(x);
    }
    integral *= h;
    return integral;
}

/* Trapezoidal rule integration */
double Trapezoidal(double local_a, double local_b, int local_n, double h, double (*f)(double x)) {
    double integral = (f(local_a) + f(local_b)) / 2.0;
    double x = local_a + h;
    for (int i = 1; i < local_n; i++, x += h) {
        integral += f(x);
    }
    integral *= h;
    return integral;
}

/* Simpson's rule integration */
double Simpsons(double local_a, int local_n, double h, double (*f)(double x)) {
    double integral = 0.0;
    double x = local_a;
    for (int i = 0; i < local_n; i++) {
        integral += f(x) + 4 * f(x + h / 2) + f(x + h);
        x += h;
    }
    integral *= h / 6;
    return integral;
}

/* Optimized Rectangle rule integration */
double Rectangle2(double local_a, int local_n, double h, double (*f)(double x)) {
    double integral = 0.0;
    double x = local_a;
    for (int i = 0; i < local_n; i++) {
        x += h;
        integral += f(x);
    }
    integral *= h;
    return integral;
}

/* Modified Trapezoidal rule integration */
double Trapezoidal2(double local_a, double local_b, int local_n, double h, double (*f)(double x)) {
    double integral = (f(local_a) + f(local_b)) / 2.0;
    double x = local_a;
    for (int i = 0; i < local_n; i++) {
        x += h;
        integral += f(x);
    }
    integral *= h;
    return integral;
}


/* Mathematical function pi(x) = 4 / (1 + x^2) */
double pi(double x) {
    return 4.0 / (1.0 + x * x);
}

/* Mathematical function ex2(x) = exp(-x^2) */
double ex2(double x) {
    return exp(-x * x);
}

/* Mathematical function sinx(x) = sin(x) / x */
double sinx(double x) {
    if (fabs(x) < Z0) { // Handling the singularity at x = 0
        return 1.0;
    }
    return sin(x) / x;
}

/* Mathematical function x2(x) = 3 * x^2 */
double x2(double x) {
    return 3.0 * x * x;
}


double select_integral_method(int rule_flag, int fun_flag, double local_a, double local_b, int local_n, double h) {
    if (rule_flag == 1 && fun_flag == 1) return Rectangle(local_a, local_n, h, pi);
    else if (rule_flag == 1 && fun_flag == 2) return Rectangle(local_a, local_n, h, ex2);
    else if (rule_flag == 1 && fun_flag == 3) return Rectangle(local_a, local_n, h, sinx);
    else if (rule_flag == 1 && fun_flag == 4) return Rectangle(local_a, local_n, h, x2);
    else if (rule_flag == 2 && fun_flag == 1) return Trapezoidal(local_a, local_b, local_n, h, pi);
    else if (rule_flag == 2 && fun_flag == 2) return Trapezoidal(local_a, local_b, local_n, h, ex2);
    else if (rule_flag == 2 && fun_flag == 3) return Trapezoidal(local_a, local_b, local_n, h, sinx);
    else if (rule_flag == 2 && fun_flag == 4) return Trapezoidal(local_a, local_b, local_n, h, x2);
    else if (rule_flag == 3 && fun_flag == 1) return Simpsons(local_a, local_n, h, pi);
    else if (rule_flag == 3 && fun_flag == 2) return Simpsons(local_a, local_n, h, ex2);
    else if (rule_flag == 3 && fun_flag == 3) return Simpsons(local_a, local_n, h, sinx);
    else if (rule_flag == 3 && fun_flag == 4) return Simpsons(local_a, local_n, h, x2);
    // Add cases for the remaining combinations...
    else return 0;  // Default to 0 if no valid combination found
}