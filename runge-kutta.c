#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define STEPSIZE 0.001
#define TIMESTEPS 1000000
#define DIMENSION 2

#define IZHIKEVICH_I 0.0
#define IZHIKEVICH_A 0.02
#define IZHIKEVICH_B 2.0
#define IZHIKEVICH_C -30.0
#define IZHIKEVICH_D 4.0

#define EVENT_PRECISION 0.00001

struct Point
{
    long double data[DIMENSION];
    long double time;
};

int log_msg(
    const char * format,
    ...
    );

int save_data(
    struct Point * steps
    );

struct Point ordinary_differential_equations(
    double time,
    struct Point point
    );

struct Point * rk_method(
    struct Point t_0
    );

int log_msg(
    const char * format,
    ...
    )
{
/**
Message Logger
log_msg("INFO : argc = %d.", argc);
log_msg("INFO : argv[0] = %s.", argv[0]);
**/
    va_list args;
    FILE * fp = NULL;
    fp = fopen("rk_log.txt", "a+");
    assert(fp);
    va_start(args, format);
    vfprintf(fp, format, args);
    va_end(args);
    fclose(fp);
    return 0;
}

int save_data(
    struct Point * step 
    )
{
/**
Save data to file
in comma separated values 
(csv) format.
**/
    FILE *fp;
    fp = NULL;
    fp = fopen("rk_out.csv", "w+");
    assert(fp);
    for(int i=0; i<DIMENSION + 1; i++){
        if(0 == i)
            fprintf(fp, "t, ");  
        else
            fprintf(fp, "d_%d, ", i-1);  
    }
    fprintf(fp, "I, ");  
    fprintf(fp, "a, ");  
    fprintf(fp, "b, ");  
    fprintf(fp, "c, ");  
    fprintf(fp, "d, ");  
    fprintf(fp, "STEPSIZE, ");  
    fprintf(fp, "TIMESTEPS, ");  
    fprintf(fp, "DIMENSION, ");  
    for(int t=0; t<TIMESTEPS; t++){
        fprintf(fp, "\n");
        for(int i=0; i<DIMENSION + 1; i++){
            if(0 == i)
                fprintf(fp, "%Lf, ", step[t].time);  
            else{
                if(isnan(step[t].data[i-1]))
                    fprintf(fp, "NaN, ");  
                else
                    fprintf(fp, "%Lf, ", step[t].data[i-1]);  
            }
        }
        fprintf(fp, "%f, ",IZHIKEVICH_I);  
        fprintf(fp, "%f, ",IZHIKEVICH_A);  
        fprintf(fp, "%f, ",IZHIKEVICH_B);  
        fprintf(fp, "%f, ",IZHIKEVICH_C);  
        fprintf(fp, "%f, ",IZHIKEVICH_D);

        fprintf(fp, "%f, ", STEPSIZE);  
        fprintf(fp, "%d, ", TIMESTEPS);  
        fprintf(fp, "%d, ", DIMENSION);  
    }
    fclose(fp);
    return 0;
}

struct Point ordinary_differential_equations(
    double time,
    struct Point point 
    )
{
/**
The ODE system to be integrated.

TEST SYSTEM

dv/dt = 2v - 3u,
du/dt = 2v - 2u.

assert(2 == DIMENSION);
rate_of_change.data[0] = 2 * point.data[0] - 3 * point.data[1];
rate_of_change.data[1] = 2 * point.data[0] - 2 * point.data[1];

Note, trace = 0, determinant > 0,
=> fixed point is neutrally stable center.


IZHIKEVICH NEURON MODEL

dv/dt = 0.04v^2 + 5v + 140 - u + I,
du/dt = a(bv - u),
if v >= 30,
then v -> c,
and u -> u + d.

Note, the system is autonomous.
**/
    struct Point rate_of_change;
    assert(0 == isnan(time));
    assert(2 == DIMENSION);
    rate_of_change.data[0] = 0.04 * point.data[0] * point.data[0] + 5 * point.data[0] + 140 - point.data[1] + IZHIKEVICH_I;
    rate_of_change.data[1] = IZHIKEVICH_A * (IZHIKEVICH_B * point.data[0] - point.data[1]);

    return rate_of_change;
}

struct Point rk_step(
    struct Point s,
    long double size,
    int t
    )
{
    struct Point this_step;
    struct Point next_step;
    struct Point k_1;
    struct Point k_2;
    struct Point k_3;
    struct Point k_4;

    //k_1 = h*f(t_n, x_n),
    for(int i = 0; i < DIMENSION; i++){
        this_step.data[i] = s.data[i];
    }
    k_1 = ordinary_differential_equations(t, this_step);
    for(int i = 0; i < DIMENSION; i++){
        k_1.data[i] = size * k_1.data[i];
    }
    //k_2 = h*f(t_n + h / 2, x_n + k_1 / 2)
    for(int i = 0; i < DIMENSION; i++){
        this_step.data[i] = s.data[i] + k_1.data[i] / 2;
    }
    k_2 = ordinary_differential_equations(t + size / 2, this_step);
    for(int i = 0; i < DIMENSION; i++){
        k_2.data[i] = size * k_2.data[i];
    }
    //k_3 = h*f(t_n + h / 2, x_n + k_2 / 2)
    for(int i = 0; i < DIMENSION; i++){
        this_step.data[i] = s.data[i] + k_2.data[i] / 2;
    }
    k_3 = ordinary_differential_equations(t + size / 2, this_step);
    for(int i = 0; i < DIMENSION; i++){
        k_3.data[i] = size * k_3.data[i];
    }
    //k_4 = h*f(t_n + h, x_n + k_3)
    for(int i = 0; i < DIMENSION; i++){
        this_step.data[i] = s.data[i] + k_3.data[i];
    }
    k_4 = ordinary_differential_equations(t + size, this_step);
    for(int i = 0; i < DIMENSION; i++){
        k_4.data[i] = size * k_4.data[i];
    }
    //x_{n+1} = x_n + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
    for(int i = 0; i < DIMENSION; i++){
        next_step.data[i] = this_step.data[i] + (k_1.data[i] + 2 * k_2.data[i] + 2 * k_3.data[i] + k_4.data[i]) / 6;
    }
    return next_step;
}
   
struct Point * rk_method(
    struct Point t_0
    )
{
/**
Runge-Kutta fourth order method
k_1 = h*f(t_n, x_n),
k_2 = h*f(t_n + h / 2, x_n + k_1 / 2)
k_3 = h*f(t_n + h / 2, x_n + k_2 / 2)
k_4 = h*f(t_n + h, x_n + k_3)
t_{n+1} = t_n + h
x_{n+1} = x_n + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
**/
    struct Point * step;
    step = NULL;
    step = malloc(TIMESTEPS * sizeof(struct Point));
    assert(step);

    for(int i = 0; i < DIMENSION; i++){
        step[0].data[i] = t_0.data[i];
    }

    struct Point this_step;
    struct Point next_step;
    long double this_step_size;
    //t_{n+1} = t_n + h
    for(int t = 0; t < TIMESTEPS-1; t++){
        //reset condition
        assert(2 == DIMENSION);
        if(30.0 < step[t].data[0])
        {
            this_step.data[0] = IZHIKEVICH_C;
            this_step.data[1] = step[t].data[1] + IZHIKEVICH_D;
        }
        else
        {
            for(int i = 0; i < DIMENSION; i++){
                this_step.data[i] = step[t].data[i];
            }
        }

        //normal step
        next_step = rk_step(this_step, STEPSIZE, t);

        //event location
        assert(2 == DIMENSION);
        if(30.0 < next_step.data[0])
        {
            this_step_size = STEPSIZE;
            for(int j = 0; j < 100; j++)
            {
                this_step_size = this_step_size / 2.0;
                next_step = rk_step(this_step, this_step_size, t);
                if(30 + EVENT_PRECISION > next_step.data[0])
                {
                    break;
                }
            }
            step[t+1].time = step[t].time + this_step_size;
        }
        else 
        {
            step[t+1].time = step[t].time + STEPSIZE;
        }
        for(int i = 0; i < DIMENSION; i++){
            step[t+1].data[i] = next_step.data[i];
        }
    }
    return step;
}

int main(
    int argc,
    char * argv[]
    )
{
    log_msg("INFO: Starting up.\n");
    clock_t start;
    start = clock();
    log_msg("INFO: argc = %d\n", argc);
    for(int i=0; i<argc; i++)
        log_msg("INFO: argv[%d] = %s\n", i, argv[i]);

    log_msg("INFO: I = %f\n", IZHIKEVICH_I);
    log_msg("INFO: a = %f\n", IZHIKEVICH_A);
    log_msg("INFO: b = %f\n", IZHIKEVICH_B);
    log_msg("INFO: c = %f\n", IZHIKEVICH_C);
    log_msg("INFO: d = %f\n", IZHIKEVICH_D);

    log_msg("INFO: STEPSIZE = %f\n", STEPSIZE);
    log_msg("INFO: TIMESTEPS = %d\n", TIMESTEPS);
    log_msg("INFO: DIMENSION = %d\n", DIMENSION);

    //initial conditions
    assert(3 == argc);
    struct Point t_0;
    long double arg;
    for(int i = 0; i < argc - 1; i++)
    {
        sscanf(argv[1 + i], "%Lf", &arg);
        t_0.data[i] = arg;
        log_msg("INFO: t_0[%d] = %Lf\n", i, t_0.data[i]);
    }
    clock_t stop;
    stop = clock();
    long double total;
    total = (long double) (stop - start) / CLOCKS_PER_SEC;
    log_msg("INFO: Done in %Lf seconds.\n", total);

    //RK4
    log_msg("INFO: Running Runge-Kutta fourth order method.\n");
    start = clock();
    struct Point * step;
    step = NULL;
    step = rk_method(t_0);
    assert(step);
    stop = clock();
    total = (long double) (stop - start) / CLOCKS_PER_SEC;
    log_msg("INFO: Done in %Lf seconds.\n", total);

    //save data
    log_msg("INFO: Saving to file.\n");
    start = clock();
    save_data(step);
    free(step);
    stop = clock();
    total = (long double) (stop - start) / CLOCKS_PER_SEC;
    log_msg("INFO: Done in %Lf seconds.\n", total);
    log_msg("INFO: End.\n");
    return 0;
}

