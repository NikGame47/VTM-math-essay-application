#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#define _USE_MATH_DEFINES

void muller_sphere_surface(int n, double *point) {
    double norm_sq = 0.0;
    
    // Шаг 1: Генерация n независимых нормальных величин
    for (int i = 0; i < n; i++) {
        // Box-Muller transform для N(0,1)
        double u1 = 1.0 - rand() / (RAND_MAX + 1.0);  // (0, 1]
        double u2 = 1.0 - rand() / (RAND_MAX + 1.0);  // (0, 1]
        
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        point[i] = z;
        norm_sq += z * z;
    }
    
    // Шаг 2: Нормировка на единичную длину
    double norm = sqrt(norm_sq);
    for (int i = 0; i < n; i++) {
        point[i] /= norm;
    }

    point[n]= 1.0;
}

void harmon_major_sphere_surface(int n, double *point) {
    // Генерируем n+1 экспоненциальных случайных величин
    double *E;
    E = malloc(sizeof(double) * (n + 1));
    double sum = 0.0;
    
    for (int i = 0; i <= n; i++) {
        // Экспоненциальное распределение: E = -ln(U), U ~ Uniform(0,1]
        E[i] = -log(1.0 - rand() / (RAND_MAX + 1.0));
        sum += E[i];
    }
    
    // Берем первые n координат: X_i = sqrt(E_i / sum)
    for (int i = 0; i < n; i++) {
        point[i] = sqrt(E[i] / sum);
    }
    
    // Добавляем случайные знаки для покрытия всей сферы
    for (int i = 0; i < n; i++) {
        if (rand() % 2) {
            point[i] = -point[i];
        }
    }
    
    point[n] = 1.0;
    
    free(E);
}

void polar_sphere_2d(double *point) {
    // Генерируем равномерный угол
    double theta = 2.0 * M_PI * rand() / (RAND_MAX + 1.0);
    
    // Координаты на единичной окружности
    point[0] = cos(theta);
    point[1] = sin(theta);
    point[2] = 1.0;
}

void polar_sphere_3d(double *point) {
    // Генерируем равномерные углы
    double u = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;  // cos(theta) ~ Uniform[-1,1]
    double phi = 2.0 * M_PI * rand() / (RAND_MAX + 1.0);  // азимутальный угол
    
    // Вычисляем sin(theta)
    double sin_theta = sqrt(1.0 - u * u);
    
    // Декартовы координаты и норма
    point[0] = sin_theta * cos(phi);
    point[1] = sin_theta * sin(phi);
    point[2] = u;
    point[3] = 1.0;
}

void quaternion_sphere_4d(double *point) {
    // Генерация трех равномерных случайных величин
    double u1 = rand() / (RAND_MAX + 1.0);  // (0, 1)
    double u2 = rand() / (RAND_MAX + 1.0);  // (0, 1)
    double u3 = rand() / (RAND_MAX + 1.0);  // (0, 1)
    
    // Параметризация равномерного кватерниона
    double r1 = sqrt(1.0 - u1);
    double r2 = sqrt(u1);
    double theta1 = 2.0 * M_PI * u2;
    double theta2 = 2.0 * M_PI * u3;
    
    // Компоненты кватерниона (точка на S3) и норма
    point[0] = r1 * sin(theta1);  // q0 = x
    point[1] = r1 * cos(theta1);  // q1 = y  
    point[2] = r2 * sin(theta2);  // q2 = z
    point[3] = r2 * cos(theta2);  // q3 = w
    point[4] = 1.0;
}

void rejection_projection_sphere(int n, double *point, int isonsphere) {
    double norm_sq;
    
    // Генерация точки внутри куба [-1,1]^n
    do {
        norm_sq = 0.0;
        for (int i = 0; i < n; i++) {
            point[i] = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
            norm_sq += point[i] * point[i];
        }
    } while (norm_sq > 1.0 || norm_sq == 0.0);

    point[n] = norm_sq;
    
    if (!isonsphere) return;

    // Проекция на единичную сферу
    double norm = sqrt(norm_sq);
    for (int i = 0; i < n; i++) {
        point[i] /= norm;
    }

    point[n] = 1.0;
}

void move_into_sphere(int n, double *point) {
    // Генерация радиуса
    double u = 1.0 - rand() / (RAND_MAX + 1.0);
    
    // Оптимизация для часто используемых размерностей
    double radius;
    switch (n) {
        case 2:
            radius = sqrt(u);  
            break;
        case 3:
            radius = cbrt(u);  
            break;
        case 4:
            radius = sqrt(sqrt(u));  
            break;
        default:
            radius = pow(u, 1.0 / n);
            break;
    }
    
    // Умножение координат на радиус
    for (int i = 0; i < n; i++) {
        point[i] *= radius;
    }

    point[n] = radius * radius; // Норма точки
}

double doyouknowdaway(int n, double *point, double *vect) {
    double D, b, c;
    b = 0.0;
    c = point[n] - 1.0;
    for (int i = 0; i < n; i++) {
        b += point[i] * vect[i]; 
    }
    D = b * b - c; // D/4 = b^2/4 - 4ac; b == 2(x, u), a == 1 => пусть b = (x, u) и D = b^2 - c
    return sqrt(D) - b;
}

//argv[1] - размерность, argv[2] - метод, argv[3] - количество точек, argv[4] - напрвление (1 - случайное)
//Мюллер = 0, Хармон-Майор = 1, полярная = 2, сферическая = 3, кватерион = 4, отбраковка = 5
int main(int argc, char **argv) {
    double ans = 0.0, *point, *vect;
    int n, method, count, israndom;

    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &method);
    sscanf(argv[3], "%d", &count);
    sscanf(argv[4], "%d", &israndom);

    point = malloc(sizeof(double) * (n + 1));
    vect = malloc(sizeof(double) * (n + 1));
    for (int i = 1; i < n; i++) {
        vect[i] = 0.0;
    }
    vect[0] = 1.0;
    vect[n] = 1.0;
    
    for (int i = 0; i < count; i++) {
        switch (method) {
            case 0:
                muller_sphere_surface(n, point);
                move_into_sphere(n, point);
                if (israndom) muller_sphere_surface(n, vect);
                break;
            case 1:
                harmon_major_sphere_surface(n, point);
                move_into_sphere(n, point);
                if (israndom) harmon_major_sphere_surface(n, vect);
                break;
            case 2:
                polar_sphere_2d(point);
                move_into_sphere(n, point);
                if (israndom) polar_sphere_2d(vect);
                break;
            case 3:
                polar_sphere_3d(point);
                move_into_sphere(n, point);
                if (israndom) polar_sphere_3d(vect);
                break;
            case 4:
                quaternion_sphere_4d(point);
                move_into_sphere(n, point);
                if (israndom) quaternion_sphere_4d(vect);
                break;
            case 5:
                rejection_projection_sphere(n, point, 0);
                if (israndom) rejection_projection_sphere(n, vect, 1);
                break;
        }

        ans += doyouknowdaway(n, point, vect);

        if (i == 10 || i == 100 || i == 1000 || i == 10000 || i == 100000 || i == 1000000) {
            printf("%d %lf\n", i, ans / (i + 1));
        }
    }

    return 0;
}
