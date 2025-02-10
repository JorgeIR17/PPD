#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>  // Para manipular pipes
#include <sys/time.h> // Para usleep
#include <omp.h> // Anadir -fopenmp al compilar
#include <mpi.h>

//mpicc -fopenmp -o ACOparalelo ACOparalelo.c -lm
//mpirun -np 4 ./ACOparalelo

#define NUM_HORMIGAS 100
#define MAX_ITERACIONES 10
#define ALPHA 1.0
#define BETA 2.0
#define EVAPORACION 0.5
#define Q 100.0

int NUM_PUNTOS = 0;
int PUNTO_INICIAL = 0;
int PUNTO_FINAL = 0;

// Funcion para inicializar la matriz de feromonas
void inicializar_feromonas(double ***feromonas) 
{
    *feromonas = (double **)malloc(NUM_PUNTOS * sizeof(double *));
    if (*feromonas == NULL) 
    {
        printf("error al reservar matriz\n");
        perror("Error al reservar memoria para la matriz de feromonas");
    }
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        (*feromonas)[i] = (double *)malloc(NUM_PUNTOS * sizeof(double));
        if ((*feromonas)[i] == NULL) 
        {
            printf("error al reservar matriz\n");
            perror("Error al reservar memoria para la matriz de feromonas");
        }

        for (int j = 0; j < NUM_PUNTOS; j++)
        {
            (*feromonas)[i][j] = 1.0; // Inicializar feromonas
        }
    }
}

// Funcion para liberar la memoria de la matriz de feromonas
void liberar_feromonas(double **feromonas) 
{
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {                
        free(feromonas[i]);
    }
    free(feromonas);
}

// Funcion para enviar un comando al servidor y recibir la respuesta
void comunicacion(char* pipe_wd, char* pipe_rd, const char *comando, const char *datos, char *respuesta) 
{
    char mensaje[256];
    snprintf(mensaje, sizeof(mensaje), "%s|%s", comando, datos);
    usleep(200); //Espera para enviar
    int write_fd = open(pipe_wd, O_WRONLY);
    if (write_fd == -1) 
    {
        perror("Error al abrir pipe_wd");
    }
    if (write(write_fd, mensaje, strlen(mensaje) + 1) == -1) 
    {
        printf("Error al escribir\n");
        perror("Error al escribir en el pipe");
    }
    close(write_fd);

    int read_fd = open(pipe_rd, O_RDONLY);
    if (read_fd == -1) 
    {
        perror("Error al abrir pipe_rd");
    }
    if (read(read_fd, respuesta, 256) == -1) 
    {
        perror("Error al leer del pipe");
        printf("Error al leer\n");
    }
    close(read_fd);
}


// Dadas las coordenadas del punto, me da el numero de punto dentro del rectangulo que define el oceano
int obtener_numeracion(float y, float x, float x_min, float x_max, float y_min, float y_max, float resolucion) 
{
    // Calcular puntos por fila
    int puntos_por_fila = (int)floor((x_max - x_min) / resolucion) + 1;

    // Determinar índices del punto
    int indice_fila = (int)floor((y - y_min) / resolucion);    // Indice de la fila
    int indice_columna = (int)floor((x - x_min) / resolucion); // Indice de la columna

    // Calcular la numeracion del punto
    int numero_punto = (indice_fila * puntos_por_fila) + indice_columna;

    return numero_punto;
}


// Dado un punto (numero), me devuelve sus coordenadas de longitud y latitud (x, y)
void obtener_coordenadas(int numero, float x_min, float x_max, float y_min, float y_max, float resolucion, float *x, float *y) 
{
    // Calcular puntos por fila
    int puntos_por_fila = (int)floor((x_max - x_min) / resolucion) + 1;

    // Calcular fila y columna del punto
    int fila = floor(numero / puntos_por_fila);
    int columna = numero % puntos_por_fila;

    // Calcular coordenadas
    *x = x_min + (columna * resolucion);
    *y = y_min + (fila * resolucion);
}

// Funcion para obtener el tiempo entre dos puntos
float tiempo(int origen, int destino, char* pipe_wd, char* pipe_rd) 
{
    float time = 0.0;
    char datos[100], respuesta[256] = "";
    snprintf(datos, sizeof(datos), "%d %d", origen, destino);
    comunicacion(pipe_wd, pipe_rd, "calcular_tiempo", datos, respuesta);

    // Verificar la respuesta
    if (strncmp(respuesta, "Error", 5) == 0) 
    {
        fprintf(stderr, "Error del servidor: %s\n", respuesta);
        return -1.0;
    }
    time = atof(respuesta);
    //printf("%.2f\n", time);
    if(isnan(time))
        return -1.0;
    return time;
}

// Funcion para calcular la longitud de una ruta
double calcular_longitud_ruta(int *ruta, char* pipe_wd, char* pipe_rd) 
{
    double longitud = 0.0;
    int i = 0;

    while (i < NUM_PUNTOS-1 && ruta[i] != -1 && ruta[i+1] != -1) // Descarta -1
    { 
        #pragma omp critical
        longitud += tiempo(ruta[i], ruta[i + 1], pipe_wd, pipe_rd);
        i++;
    }
    
    return longitud;
}

// Funcion para seleccionar la siguiente ciudad
int seleccionar_ciudad(double **feromonas, int ciudad_actual, int *visitado, char* pipe_wd, char* pipe_rd) 
{
    double *probabilidades = (double *)calloc(NUM_PUNTOS, sizeof(double));
    double suma = 0.0;
    int ultimo_no_visitado = -1;
    float t = 0.0;
    double heuristica = 0.0;

    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        if (!visitado[i]) 
        {
            #pragma omp critical
            t = tiempo(ciudad_actual, i, pipe_wd, pipe_rd);
            heuristica = 1.0 / (t + 1e-6);
            probabilidades[i] = pow(feromonas[ciudad_actual][i], ALPHA) * pow(heuristica, BETA);
            suma += probabilidades[i];
            ultimo_no_visitado = i;
        } 
        else 
            probabilidades[i] = 0.0;
    }

    // Normalizar probabilidades
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        probabilidades[i] /= suma;
        if (i > 0) 
            probabilidades[i] += probabilidades[i - 1];
    }
    probabilidades[ultimo_no_visitado] = 1.0;

    // Seleccionar ciudad
    double r = 0.0;
    r = rand() / (double)RAND_MAX;
    for (int i = 0; i < NUM_PUNTOS; i++)  
    {
        if (r <= probabilidades[i]) 
        {
            free(probabilidades);
            return i;
        }
    }

    free(probabilidades);
    return -1;
}

// Funcion para actualizar las feromonas
void actualizar_feromonas(double **feromonas, int **rutas, double longitudes[], int hormigas) 
{
    #pragma omp parallel for collapse(2) shared(feromonas)
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        for (int j = 0; j < NUM_PUNTOS; j++) 
        {
            feromonas[i][j] *= (1 - EVAPORACION);
        }
    }

    #pragma omp parallel for shared(feromonas, rutas, longitudes)
    for (int k = 0; k < hormigas; k++) 
    {
        double contribucion = Q / longitudes[k];
        for (int i = 0; i < NUM_PUNTOS - 1 && rutas[k][i] != -1 && rutas[k][i + 1] != -1; i++) 
        {
            int c1 = rutas[k][i], c2 = rutas[k][i + 1];

            #pragma omp atomic
            feromonas[c1][c2] += contribucion;

            #pragma omp atomic
            feromonas[c2][c1] += contribucion;
        }
    }
}


// Función principal
int main(int argc, char** argv) 
{
    int rank, size;

    int n_hormigas1 = 0;
    int n_hormigas2 = 0;
    int n_hormigas3 = 0;
    int n_hormigas4 = 0;
    int n_hormigas5 = 0;
    int n_hormigas6 = 0;
    int n_hormigas7 = 0;
    int n_hormigas8 = 0;
    int n_hormigas9 = 0;
    int n_hormigas10 = 0;
    int n_hormigas11 = 0;
    int n_hormigas12 = 0;
    int n_hormigas13 = 0;
    int n_hormigas14 = 0;
    int n_hormigas15 = 0;
    int n_hormigas16 = 0;

    FILE *f = fopen("nodes.txt", "r");
    if(f == NULL)
    {
        perror("Error al abrir el archivo\n");
    }

    if(fscanf(f, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &n_hormigas1, &n_hormigas2, &n_hormigas3, &n_hormigas4, &n_hormigas5, &n_hormigas6, 
    &n_hormigas7, &n_hormigas8, &n_hormigas9, &n_hormigas10, &n_hormigas11, &n_hormigas12, &n_hormigas13, &n_hormigas14, &n_hormigas15, &n_hormigas16) != 16)
    {
        fprintf(stderr, "Error al leer los datos del archivo\n");
        fclose(f);
        return 1;
    }

    if(n_hormigas1 + n_hormigas2 + n_hormigas3 + n_hormigas4 + n_hormigas5 + n_hormigas6 + 
    n_hormigas7 + n_hormigas8 + n_hormigas9 + n_hormigas10 + n_hormigas11 + n_hormigas12 + n_hormigas13 + n_hormigas14 + n_hormigas15 + n_hormigas16 != 100)
    {
        fprintf(stderr, "Los datos son incorrectos\n");
        fclose(f);
        return 1;
    }

    n_hormigas1 = (NUM_HORMIGAS * n_hormigas1) / 100;
    n_hormigas2 = (NUM_HORMIGAS * n_hormigas2) / 100;
    n_hormigas3 = (NUM_HORMIGAS * n_hormigas3) / 100;
    n_hormigas4 = (NUM_HORMIGAS * n_hormigas4) / 100;
    n_hormigas5 = (NUM_HORMIGAS * n_hormigas5) / 100;
    n_hormigas6 = (NUM_HORMIGAS * n_hormigas6) / 100;
    n_hormigas7 = (NUM_HORMIGAS * n_hormigas7) / 100;
    n_hormigas8 = (NUM_HORMIGAS * n_hormigas8) / 100;
    n_hormigas9 = (NUM_HORMIGAS * n_hormigas9) / 100;
    n_hormigas10 = (NUM_HORMIGAS * n_hormigas10) / 100;
    n_hormigas11 = (NUM_HORMIGAS * n_hormigas11) / 100;
    n_hormigas12 = (NUM_HORMIGAS * n_hormigas12) / 100;
    n_hormigas13 = (NUM_HORMIGAS * n_hormigas13) / 100;
    n_hormigas14 = (NUM_HORMIGAS * n_hormigas14) / 100;
    n_hormigas15 = (NUM_HORMIGAS * n_hormigas15) / 100;
    n_hormigas16 = (NUM_HORMIGAS * n_hormigas16) / 100;

    fclose(f);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_threads = 1;
    omp_set_num_threads(n_threads);

    double **feromonas = NULL;
    double *vector_feromonas = (double*)malloc(NUM_PUNTOS * NUM_PUNTOS * sizeof(double));
    double *vector_total = (double*)malloc(NUM_PUNTOS * NUM_PUNTOS * sizeof(double));

    float lat_i, lon_i, lat_f, lon_f; // Latitudes y longitudes del punto inicial y del punto final en mi trayectoria
    lat_i = -6;
    lon_i = 1;
    lat_f = 4;
    lon_f = 9;
    float lat_i_mapa, lon_i_mapa, lat_f_mapa, lon_f_mapa, resolucion_mapa; // Latitudes y longitudes que delimitan el oceano

    // Abrir pipes para comunicacion
    char* pipe_wd = "/tmp/pipe_wd";
    char* pipe_rd = "/tmp/pipe_rd";
    char* pipe_wd2 = "/tmp/pipe_wd2";
    char* pipe_rd2 = "/tmp/pipe_rd2";
    char* pipe_wd3 = "/tmp/pipe_wd3";
    char* pipe_rd3 = "/tmp/pipe_rd3";
    char* pipe_wd4 = "/tmp/pipe_wd4";
    char* pipe_rd4 = "/tmp/pipe_rd4";
    char* pipe_wd5 = "/tmp/pipe_wd5";
    char* pipe_rd5 = "/tmp/pipe_rd5";
    char* pipe_wd6 = "/tmp/pipe_wd6";
    char* pipe_rd6 = "/tmp/pipe_rd6";
    char* pipe_wd7 = "/tmp/pipe_wd7";
    char* pipe_rd7 = "/tmp/pipe_rd7";
    char* pipe_wd8 = "/tmp/pipe_wd8";
    char* pipe_rd8 = "/tmp/pipe_rd8";
    char* pipe_wd9 = "/tmp/pipe_wd9";
    char* pipe_rd9 = "/tmp/pipe_rd9";
    char* pipe_wd10 = "/tmp/pipe_wd10";
    char* pipe_rd10 = "/tmp/pipe_rd10";
    char* pipe_wd11 = "/tmp/pipe_wd11";
    char* pipe_rd11 = "/tmp/pipe_rd11";
    char* pipe_wd12 = "/tmp/pipe_wd12";
    char* pipe_rd12 = "/tmp/pipe_rd12";
    char* pipe_wd13 = "/tmp/pipe_wd13";
    char* pipe_rd13 = "/tmp/pipe_rd13";
    char* pipe_wd14 = "/tmp/pipe_wd14";
    char* pipe_rd14 = "/tmp/pipe_rd14";
    char* pipe_wd15 = "/tmp/pipe_wd15";
    char* pipe_rd15 = "/tmp/pipe_rd15";
    char* pipe_wd16 = "/tmp/pipe_wd16";
    char* pipe_rd16 = "/tmp/pipe_rd16";

    srand(time(NULL));

    char respuesta[256] = "";

    if (rank == 0)
    {
        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", n_hormigas1, n_hormigas2, n_hormigas3, n_hormigas4, n_hormigas5, n_hormigas6, 
    n_hormigas7, n_hormigas8, n_hormigas9, n_hormigas10, n_hormigas11, n_hormigas12, n_hormigas13, n_hormigas14, n_hormigas15, n_hormigas16);

        comunicacion(pipe_wd, pipe_rd, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);
        printf("Numero de puntos del mapa: %d\n", NUM_PUNTOS);

        comunicacion(pipe_wd, pipe_rd, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);
        printf("Limites del mapa recibidos del servidor: %.2f, %.2f, %.2f, %.2f, %.2f\n",
               lat_i_mapa, lon_i_mapa, lat_f_mapa, lon_f_mapa, resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);

        printf("Punto inicial y punto final de la trayectoria: %d, %d\n", PUNTO_INICIAL, PUNTO_FINAL);
    }

    if (rank == 1)
    {
        comunicacion(pipe_wd2, pipe_rd2, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd2, pipe_rd2, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 2)
    {
        comunicacion(pipe_wd3, pipe_rd3, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd3, pipe_rd3, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 3)
    {
        comunicacion(pipe_wd4, pipe_rd4, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);
        printf("Numero de puntos del mapa: %d\n", NUM_PUNTOS);

        comunicacion(pipe_wd4, pipe_rd4, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 4)
    {
        comunicacion(pipe_wd5, pipe_rd5, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);
        printf("Numero de puntos del mapa: %d\n", NUM_PUNTOS);

        comunicacion(pipe_wd5, pipe_rd5, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 5)
    {
        comunicacion(pipe_wd6, pipe_rd6, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);
        printf("Numero de puntos del mapa: %d\n", NUM_PUNTOS);

        comunicacion(pipe_wd6, pipe_rd6, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 6)
    {
        comunicacion(pipe_wd7, pipe_rd7, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd7, pipe_rd7, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 7)
    {
        comunicacion(pipe_wd8, pipe_rd8, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd8, pipe_rd8, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 8)
    {
        comunicacion(pipe_wd9, pipe_rd9, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd9, pipe_rd9, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);

    }

    if (rank == 9)
    {
        comunicacion(pipe_wd10, pipe_rd10, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd10, pipe_rd10, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 10)
    {
        comunicacion(pipe_wd11, pipe_rd11, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd11, pipe_rd11, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 11)
    {
        comunicacion(pipe_wd12, pipe_rd12, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd12, pipe_rd12, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 12)
    {
        comunicacion(pipe_wd13, pipe_rd13, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd13, pipe_rd13, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 13)
    {
        comunicacion(pipe_wd14, pipe_rd14, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd14, pipe_rd14, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 14)
    {
        comunicacion(pipe_wd15, pipe_rd15, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd15, pipe_rd15, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    if (rank == 15)
    {
        comunicacion(pipe_wd16, pipe_rd16, "num_puntos", "", respuesta);
        NUM_PUNTOS = atoi(respuesta);

        comunicacion(pipe_wd16, pipe_rd16, "limites", "", respuesta);
        sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);

        PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
        PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    }

    inicializar_feromonas(&feromonas);

    int *mejor_ruta = (int*)malloc(NUM_PUNTOS * sizeof(int));
    

    double mejor_longitud = INFINITY;


    int *mejor_ruta2 = (int*)malloc(NUM_PUNTOS * sizeof(int));
    double mejor_longitud2 = INFINITY;

    for (int iter = 0; iter < MAX_ITERACIONES; iter++)
    {
        if(rank == 0) // PROCESO 0
        {
            int **rutas = (int **)malloc(n_hormigas1 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes[n_hormigas1];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas1; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd, pipe_rd);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes[k] = calcular_longitud_ruta(rutas[k], pipe_wd, pipe_rd);

                #pragma omp critical
                if (longitudes[k] < mejor_longitud && longitudes[k] > 0)
                {
                    mejor_longitud = longitudes[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes, n_hormigas1);

            for (int k = 0; k < n_hormigas1; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 1)
        {
            int **rutas = (int **)malloc(n_hormigas2 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes2[n_hormigas2];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas2; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd2, pipe_rd2);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes2[k] = calcular_longitud_ruta(rutas[k], pipe_wd2, pipe_rd2);

                #pragma omp critical
                if (longitudes2[k] < mejor_longitud && longitudes2[k] > 0)
                {
                    mejor_longitud = longitudes2[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes2, n_hormigas2);

            for (int k = 0; k < n_hormigas2; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 2)
        {
            int **rutas = (int **)malloc(n_hormigas3 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes3[n_hormigas3];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas3; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd3, pipe_rd3);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes3[k] = calcular_longitud_ruta(rutas[k], pipe_wd3, pipe_rd3);

                #pragma omp critical
                if (longitudes3[k] < mejor_longitud && longitudes3[k] > 0)
                {
                    mejor_longitud = longitudes3[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes3, n_hormigas3);

            for (int k = 0; k < n_hormigas3; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 3)
        {
            int **rutas = (int **)malloc(n_hormigas4 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes4[n_hormigas4];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas4; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd4, pipe_rd4);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes4[k] = calcular_longitud_ruta(rutas[k], pipe_wd4, pipe_rd4);

                #pragma omp critical
                if (longitudes4[k] < mejor_longitud && longitudes4[k] > 0)
                {
                    mejor_longitud = longitudes4[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes4, n_hormigas4);

            for (int k = 0; k < n_hormigas4; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 4)
        {
            int **rutas = (int **)malloc(n_hormigas5 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes5[n_hormigas5];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas5; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd5, pipe_rd5);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes5[k] = calcular_longitud_ruta(rutas[k], pipe_wd5, pipe_rd5);

                #pragma omp critical
                if (longitudes5[k] < mejor_longitud && longitudes5[k] > 0)
                {
                    mejor_longitud = longitudes5[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes5, n_hormigas5);

            for (int k = 0; k < n_hormigas5; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 5)
        {
            int **rutas = (int **)malloc(n_hormigas6 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes6[n_hormigas6];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas6; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd6, pipe_rd6);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes6[k] = calcular_longitud_ruta(rutas[k], pipe_wd6, pipe_rd6);

                #pragma omp critical
                if (longitudes6[k] < mejor_longitud && longitudes6[k] > 0)
                {
                    mejor_longitud = longitudes6[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes6, n_hormigas6);

            for (int k = 0; k < n_hormigas6; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 6)
        {
            int **rutas = (int **)malloc(n_hormigas7 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes7[n_hormigas7];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas7; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd7, pipe_rd7);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes7[k] = calcular_longitud_ruta(rutas[k], pipe_wd7, pipe_rd7);

                #pragma omp critical
                if (longitudes7[k] < mejor_longitud && longitudes7[k] > 0)
                {
                    mejor_longitud = longitudes7[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes7, n_hormigas7);

            for (int k = 0; k < n_hormigas7; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 7)
        {
            int **rutas = (int **)malloc(n_hormigas8 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes8[n_hormigas8];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas8; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd8, pipe_rd8);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes8[k] = calcular_longitud_ruta(rutas[k], pipe_wd8, pipe_rd8);

                #pragma omp critical
                if (longitudes8[k] < mejor_longitud && longitudes8[k] > 0)
                {
                    mejor_longitud = longitudes8[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes8, n_hormigas8);

            for (int k = 0; k < n_hormigas8; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 8)
        {
            int **rutas = (int **)malloc(n_hormigas9 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes9[n_hormigas9];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas9; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd9, pipe_rd9);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes9[k] = calcular_longitud_ruta(rutas[k], pipe_wd9, pipe_rd9);

                #pragma omp critical
                if (longitudes9[k] < mejor_longitud && longitudes9[k] > 0)
                {
                    mejor_longitud = longitudes9[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes9, n_hormigas9);

            for (int k = 0; k < n_hormigas9; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 9)
        {
            int **rutas = (int **)malloc(n_hormigas10 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes10[n_hormigas10];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas10; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd10, pipe_rd10);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes10[k] = calcular_longitud_ruta(rutas[k], pipe_wd10, pipe_rd10);

                #pragma omp critical
                if (longitudes10[k] < mejor_longitud && longitudes10[k] > 0)
                {
                    mejor_longitud = longitudes10[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes10, n_hormigas10);

            for (int k = 0; k < n_hormigas10; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 10)
        {
            int **rutas = (int **)malloc(n_hormigas11 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes11[n_hormigas11];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas11; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd11, pipe_rd11);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes11[k] = calcular_longitud_ruta(rutas[k], pipe_wd11, pipe_rd11);

                #pragma omp critical
                if (longitudes11[k] < mejor_longitud && longitudes11[k] > 0)
                {
                    mejor_longitud = longitudes11[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes11, n_hormigas11);

            for (int k = 0; k < n_hormigas11; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 11)
        {
            int **rutas = (int **)malloc(n_hormigas12 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes12[n_hormigas12];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas12; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd12, pipe_rd12);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes12[k] = calcular_longitud_ruta(rutas[k], pipe_wd12, pipe_rd12);

                #pragma omp critical
                if (longitudes12[k] < mejor_longitud && longitudes12[k] > 0)
                {
                    mejor_longitud = longitudes12[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes12, n_hormigas12);

            for (int k = 0; k < n_hormigas12; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 12)
        {
            int **rutas = (int **)malloc(n_hormigas13 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes13[n_hormigas13];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas13; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd13, pipe_rd13);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes13[k] = calcular_longitud_ruta(rutas[k], pipe_wd13, pipe_rd13);

                #pragma omp critical
                if (longitudes13[k] < mejor_longitud && longitudes13[k] > 0)
                {
                    mejor_longitud = longitudes13[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes13, n_hormigas13);

            for (int k = 0; k < n_hormigas13; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 13)
        {
            int **rutas = (int **)malloc(n_hormigas14 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes14[n_hormigas14];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas14; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd14, pipe_rd14);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes14[k] = calcular_longitud_ruta(rutas[k], pipe_wd14, pipe_rd14);

                #pragma omp critical
                if (longitudes14[k] < mejor_longitud && longitudes14[k] > 0)
                {
                    mejor_longitud = longitudes14[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes14, n_hormigas14);

            for (int k = 0; k < n_hormigas14; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }
        
        if(rank == 14)
        {
            int **rutas = (int **)malloc(n_hormigas15 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes15[n_hormigas15];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas15; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd15, pipe_rd15);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes15[k] = calcular_longitud_ruta(rutas[k], pipe_wd15, pipe_rd15);

                #pragma omp critical
                if (longitudes15[k] < mejor_longitud && longitudes15[k] > 0)
                {
                    mejor_longitud = longitudes15[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes15, n_hormigas15);

            for (int k = 0; k < n_hormigas15; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(rank == 15)
        {
            int **rutas = (int **)malloc(n_hormigas16 * sizeof(int *));
            if (rutas == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas\n");
                exit(EXIT_FAILURE);
            }

            double longitudes16[n_hormigas16];

            #pragma omp parallel for shared(rutas, feromonas, mejor_ruta, mejor_longitud)
            for (int k = 0; k < n_hormigas16; k++)
            {
                rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
                if (rutas[k] == NULL) 
                {
                    fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                    exit(EXIT_FAILURE);
                }

                memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

                int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
                rutas[k][0] = PUNTO_INICIAL;
                visitado[PUNTO_INICIAL] = 1;

                for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
                {
                    rutas[k][i] = seleccionar_ciudad(feromonas, rutas[k][i - 1], visitado, pipe_wd16, pipe_rd16);
                    visitado[rutas[k][i]] = 1;
                }
                longitudes16[k] = calcular_longitud_ruta(rutas[k], pipe_wd16, pipe_rd16);

                #pragma omp critical
                if (longitudes16[k] < mejor_longitud && longitudes16[k] > 0)
                {
                    mejor_longitud = longitudes16[k];
                    memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
                }
                free(visitado);
            }

            actualizar_feromonas(feromonas, rutas, longitudes16, n_hormigas16);

            for (int k = 0; k < n_hormigas16; k++) 
            {
                free(rutas[k]);
            }
            free(rutas);
        }

        if(size > 1)
        {
            if((iter+1 % 2) == 0)
            {
                if(rank != 0)
                {
                    #pragma omp parallel for collapse(2) shared(feromonas, vector_feromonas)
                    for (int i = 0; i < NUM_PUNTOS; i++) 
                    {
                            for (int j = 0; j < NUM_PUNTOS; j++) 
                            {
                                int k = i * NUM_PUNTOS + j; // Índice en el vector
                                vector_feromonas[k] = feromonas[i][j];
                            }
                    }
                }

                if (rank == 0) 
                {
                    MPI_Reduce(vector_feromonas, vector_total, NUM_PUNTOS * NUM_PUNTOS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                    MPI_Bcast(vector_total, NUM_PUNTOS * NUM_PUNTOS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                }    
                    
                #pragma omp parallel for shared(feromonas, vector_total)
                for (int k = 0; k < NUM_PUNTOS * NUM_PUNTOS; k++) 
                {
                    int i = k / NUM_PUNTOS; 
                    int j = k % NUM_PUNTOS;
                    feromonas[i][j] = vector_total[k] / size;
                }
                
            }
        }
        //if(rank == 0)
        printf("Iteración %d proceso %d: Mejor longitud = %.2f\n", iter + 1, rank, mejor_longitud);
    }

    if(size > 1)
    {
        if(rank == 0)
        {
            for(int i = 1; i < size; i++)
            {
                MPI_Recv(mejor_ruta2, NUM_PUNTOS, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&mejor_longitud2, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if(mejor_longitud2 < mejor_longitud)
                {
                    mejor_longitud = mejor_longitud2;
                    memcpy(mejor_ruta, mejor_ruta2, NUM_PUNTOS * sizeof(int));
                }
            }
        }
        else
        {
            MPI_Send(mejor_ruta, NUM_PUNTOS, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&mejor_longitud, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    if(rank == 0)
    {
        printf("Mejor ruta: \n");

        float lat, lon;
        for (int i = 0; i < NUM_PUNTOS && mejor_ruta[i] != -1; i++)
        {
            printf("%d ", mejor_ruta[i]);
            obtener_coordenadas(mejor_ruta[i], lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa, &lon, &lat);
            printf("Latitud: %f, Longitud: %f\n", lat, lon);
        }
        
        printf("\nTiempo de la mejor ruta: %.2f\n", mejor_longitud);
    }

    liberar_feromonas(feromonas);
    free(vector_feromonas);
    free(vector_total);

    free(mejor_ruta);

    if(rank == 0)
        comunicacion(pipe_wd, pipe_rd, "salir", "", respuesta);
    if(rank == 1)
        comunicacion(pipe_wd2, pipe_rd2, "salir", "", respuesta);
    if(rank == 2)
        comunicacion(pipe_wd3, pipe_rd3, "salir", "", respuesta);
    if(rank == 3)
        comunicacion(pipe_wd4, pipe_rd4, "salir", "", respuesta);
    if(rank == 4)
        comunicacion(pipe_wd5, pipe_rd5, "salir", "", respuesta);
    if(rank == 5)
        comunicacion(pipe_wd6, pipe_rd6, "salir", "", respuesta);
    if(rank == 6)
        comunicacion(pipe_wd7, pipe_rd7, "salir", "", respuesta);
    if(rank == 7)
        comunicacion(pipe_wd8, pipe_rd8, "salir", "", respuesta);
    if(rank == 8)
        comunicacion(pipe_wd9, pipe_rd9, "salir", "", respuesta);
    if(rank == 9)
        comunicacion(pipe_wd10, pipe_rd10, "salir", "", respuesta);
    if(rank == 10)
        comunicacion(pipe_wd11, pipe_rd11, "salir", "", respuesta);
    if(rank == 11)
        comunicacion(pipe_wd12, pipe_rd12, "salir", "", respuesta);
    if(rank == 12)
        comunicacion(pipe_wd13, pipe_rd13, "salir", "", respuesta);
    if(rank == 13)
        comunicacion(pipe_wd14, pipe_rd14, "salir", "", respuesta);
    if(rank == 14)
        comunicacion(pipe_wd15, pipe_rd15, "salir", "", respuesta);
    if(rank == 15)
        comunicacion(pipe_wd16, pipe_rd16, "salir", "", respuesta);

    MPI_Finalize();

    return 0;
}