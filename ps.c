#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define NARGS 4
#define Z0 0.000000001

int main(int argc, char *argv[]) {
    double integral;
    double a, b;
    double h;
    int n;

    double pi(double x);
    double ex2(double x);
    double sinx(double x);
    double x2(double x);

    double Rectangle(double a, int n, double h, double (*f)(double x));
    double Trapezoidal(double a, double b, int n, double h, double (*f)(double x));
    double Simpsons(double a, int n, double h, double (*f)(double x));
    double Rectangle2(double a, int n, double h, double (*f)(double x));
    double Trapezoidal2(double a, int n, double h, double (*f)(double x));

    int ch;
    int fun_flag = 1; 
    int rule_flag = 1; 

    clock_t start, end;
    double cpu_time_used;

    time_t t_start, t_end;
    double t_cpu_time_used;

    if (argc != NARGS) {
        rule_flag = 1;
        fun_flag = 1;
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
        printf("N = number of segments (rectangles, trapezoids, or arcs)\n");
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
                default:
                    printf("Usage: %s [-rule -function]\n", argv[0]);
                    return 2;
            }
        }
    }

    a = Z0;
    b = 1;

    start = clock();
    t_start = time(NULL);
    h = (b - a) / n;

    switch (rule_flag) {
        case 1:
            integral = Rectangle(a, n, h, fun_flag == 1 ? pi : (fun_flag == 2 ? ex2 : (fun_flag == 3 ? sinx : x2)));
            break;
        case 2:
            integral = Trapezoidal(a, b, n, h, fun_flag == 1 ? pi : (fun_flag == 2 ? ex2 : (fun_flag == 3 ? sinx : x2)));
            break;
        case 3:
            integral = Simpsons(a, n, h, fun_flag == 1 ? pi : (fun_flag == 2 ? ex2 : (fun_flag == 3 ? sinx : x2)));
            break;
        case 4:
            integral = Rectangle2(a, n, h, fun_flag == 1 ? pi : (fun_flag == 2 ? ex2 : (fun_flag == 3 ? sinx : x2)));
            break;
        case 5:
            integral = Trapezoidal2(a, n, h, fun_flag == 1 ? pi : (fun_flag == 2 ? ex2 : (fun_flag == 3 ? sinx : x2)));
            break;
        default:
            integral = 0;
    }

    end = clock();
    t_end = time(NULL);

    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    t_cpu_time_used = ((double)(t_end - t_start));

    printf("\nFunctions: Rules:\n");
    printf("********* ******\n");
    printf(" f[1] = 4/((x*x)+1) r[1] = Rectangle\n");
    printf(" f[2] = e^-(x^2) r[2] = Trapezoidal\n");
    printf(" f[3] = sin(x)/x r[3] = Simpsons\n");
    printf(" f[4] = 3*x^2 r[4] = Optimized rectangle\n");
    printf(" r[5] = Modified trapezoidal\n");
    printf("\nIntegrating f[%d], using rule r[%d], from %f to %f. %d segments.\n", fun_flag, rule_flag, a, b, n);
    printf("\nIntegration result (area) = %.25f\n", integral);
    printf("\nCPU-time using clock(): %f\n", cpu_time_used);
    printf("\nCPU-time using time(): %f\n", t_cpu_time_used);

    printf("==========================================================\n");
    return 0;
}

double Rectangle(double a, int n, double h, double (*f)(double x)) {
    double integral = 0;
    for (int i = 0; i < n; i++) {
        integral += h * f(a + i * h);
    }
    return integral;
}

double Trapezoidal(double a, double b, int n, double h, double (*f)(double x)) {
    double x = a;
    double integral = (f(a) + f(b)) / 2.0;
    for (int i = 1; i <= n - 1; i++) {
        x += h;
        integral += f(x);
    }
    integral *= h;
    return integral;
}

double Simpsons(double a, int n, double h, double (*f)(double x)) {
    double integral = 0;
    for (int i = 0; i < n; i++) {
        integral += (f(a + i * h) + 4 * f(a + (i + 0.5) * h) + f(a + (i + 1) * h)) / 6 * h;
    }
    return integral;
}

double Rectangle2(double a, int n, double h, double (*f)(double x)) {
    double integral = 0;
    double x = a;
    for (int i = 0; i <= n - 1; i++) {
        x += h;
        integral += f(x);
    }
    integral *= h;
    return integral;
}

double Trapezoidal2(double a, int n, double h, double (*f)(double x)) {
    double integral = 0;
    for (int i = 0; i <= n - 1; i++) {
        integral += h * (f(a + i * h) + f(a + (i + 1) * h)) / 2;
    }
    return integral;
}

double pi(double x) {
    return 4 / ((x * x) + 1);
}

double ex2(double x) {
    return pow(M_E, -(pow(x, 2)));
}

double sinx(double x) {
    return sin(x) / x;
}

double x2(double x) {
    return 3 * (x * x);
}
