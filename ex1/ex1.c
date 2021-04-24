#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define BUFFERSIZE 1000
#define BASE_ARR_SIZE 100
#define ARR_SIZE_MULTIPLY 2
#define COMMA ','
#define MALLOC_FOR_STRING(str) (char *)(malloc((strlen(str) + 1) * sizeof(char)))
#define MALLOC_ARR(type, size) (type*)malloc((size) * sizeof(type))

int d = 0;
int k = 0;
double minus_one = -1;
int max_iter = 200;
int lines_count;
void read_text_to_lines(char ***lines);
void split_line_to_strings_nums(char *str, char ***);
void count_commas_in_str(char *str, int *d);
void turn_lines_to_vectors(int *d, int *lines_count, char ***lines, double ***vers);
int is_positive_int(const char *str);
void check_command_line_args(int *argc, char const *argv[]);
double distacne(double **a, double **b);
int closest_cluster_index(double *vector, double **centers);
void multiply_scalar(double *vector, double *scalar, double *result);
void add_vectors(double *a, double *b, double *result);
void print_final(double **mat);
void malloc_for_mat_and_set_zeros(double ***mat, int lines, int cols);
void k_means(double ***centers, double ***vectors);

int main(int argc, char const *argv[])
{
    char **lines;
    double **vectors;
    double **centers;
    check_command_line_args(&argc, argv);

    read_text_to_lines(&lines);

    turn_lines_to_vectors(&d, &lines_count, &lines, &vectors);
    free(lines);

    malloc_for_mat_and_set_zeros(&centers, k, d);
    k_means(&centers, &vectors);

    print_final(centers);

    return 0;
}

void read_text_to_lines(char ***lines)
{
    int arrlen;
    char buffer[BUFFERSIZE];
    lines_count = 0;
    arrlen = BASE_ARR_SIZE;
    *lines = MALLOC_ARR(char *,arrlen);
    assert((*lines) != NULL);
    while (fgets(buffer, BUFFERSIZE, stdin) != NULL)
    {
        if (lines_count == arrlen)
        {
            arrlen *= ARR_SIZE_MULTIPLY;
            *lines = (char **)(realloc(*lines, arrlen * sizeof(char *)));
            assert(*lines != NULL);
        }
        if (d == 0)
        {
            count_commas_in_str(buffer, &d);
        }

        (*lines)[lines_count] = (char *)(malloc((strlen(buffer) + 1) * sizeof(char)));
        assert((*lines)[lines_count] != NULL);
        strcpy((*lines)[lines_count], buffer);
        (lines_count)++;
    }
}

void count_commas_in_str(char *str, int *d)
{
    int count = 0;
    size_t i;
    for (i = 0; i < strlen(str); i++)
    {
        if (str[i] == COMMA)
        {
            count++;
        }
    }
    *d = count + 1;
}

void split_line_to_strings_nums(char *str, char ***result)
{
    char *num_str;
    int j = 0;
    *result = (char **)malloc(d * sizeof(char *));
    num_str = strtok(str, ",");
    while (num_str)
    {
        (*result)[j] = MALLOC_FOR_STRING(num_str);
        strcpy((*result)[j], num_str);
        j++;
        num_str = strtok(NULL, ",");
    }
}

void turn_lines_to_vectors(int *d, int *lines_count, char ***lines, double ***vers)
{
    int i, j;
    char **data;
    char **current_line;
    (*vers) = (double **)(malloc((*lines_count) * sizeof(double *)));
    data = (char **)(malloc)((*lines_count) * sizeof(char *));

    for (i = 0; i < *lines_count; i++)
    {
        (*vers)[i] = (double *)(malloc((*d) * sizeof(double)));
        data[i] = (*lines)[i];
        data[i][strcspn(data[i], "\n")] = '\0';
        split_line_to_strings_nums(data[i], &current_line);
        for (j = 0; j < *d; j++)
        {
            char *eptr;
            (*vers)[i][j] = strtod(current_line[j], &eptr);
        }
    }
}

int is_positive_int(const char *str)
{
    char digit = str[0];
    size_t i;
    if (digit <= '0' || digit > '9')
    {
        return 0;
    }
    for (i = 1; i < strlen(str); i++)
    {
        digit = str[i];
        if (digit < '0' || digit > '9')
        {
            return 0;
        }
    }
    return 1;
}

void check_command_line_args(int *argc, char const *argv[])
{
    if (*argc == 1)
    {
        printf("k is not provided \n");
        assert(*argc > 1);
    }
    if (!is_positive_int(argv[1]))
    {
        printf("k value is not valid \n");
        assert(is_positive_int(argv[1]));
    }
    k = atoi(argv[1]);
    if (*argc > 3)
    {
        printf("Too many arguments \n");
        assert(*argc <= 2);
    }
    if (*argc == 3)
    {
        if (!is_positive_int(argv[2]))
        {
            printf("max_iter value is not valid \n");
            assert(is_positive_int(argv[2]));
        }
        else
        {
            max_iter = atoi(argv[2]);
        }
    }
}

double distacne(double **a, double **b)
{
    double sum = 0;
    int i;
    for (i = 0; i < d; i++)
    {
        sum += ((*a)[i] - (*b)[i]) * ((*a)[i] - (*b)[i]);
    }
    return sum;
}

int closest_cluster_index(double *vector, double **centers)
{
    int min_dist_index = 0, i;
    double min_dist;
    min_dist = distacne(&vector, &centers[0]);
    for (i = 0; i < k; i++)
    {
        double current_dist;
        current_dist = distacne(&vector, &centers[i]);
        if (current_dist < min_dist)
        {
            min_dist = current_dist;
            min_dist_index = i;
        }
    }
    return min_dist_index;
}

void multiply_scalar(double *vector, double *scalar, double *result)
{
    int i;
    for (i = 0; i < d; i++)
    {
        result[i] = vector[i] * (*scalar);
    }
}

void add_vectors(double *a, double *b, double *result)
{
    int i;
    for (i = 0; i < d; i++)
    {
        result[i] = a[i] + b[i];
    }
}

void print_final(double **mat)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d - 1; j++)
        {
            printf("%.4f,", mat[i][j]);
        }
        printf("%.4f", mat[i][d - 1]);
        printf("\n");
    }
}

void malloc_for_mat_and_set_zeros(double ***mat, int lines, int cols)
{
    int i;
    *mat = MALLOC_ARR(double *, lines);
    for (i = 0; i < lines; i++)
    {
        (*mat)[i] = (double *)calloc(cols, sizeof(double));
    }
}

void k_means(double ***centers, double ***vectors)
{
    int i, j, l;
    double **cluster_sums;
    double *negative_vec;
    int *cluster_size;
    int *which_cluster;
    negative_vec = MALLOC_ARR(double, d);
    malloc_for_mat_and_set_zeros(&cluster_sums, k, d);
    cluster_size = MALLOC_ARR(int, k);
    which_cluster = MALLOC_ARR(int, lines_count);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            (*centers)[i][j] = (*vectors)[i][j];
        }
        cluster_size[i] = 0;
    }
    for (i = 0; i < lines_count; i++)
    {
        which_cluster[i] = -1;
    }
    for (i = 0; i < max_iter; i++)
    {
        int differ = 0;
        for (j = 0; j < lines_count; j++)
        {
            int new_closest_cluster;
            new_closest_cluster = closest_cluster_index((*vectors)[j], *centers);
            if (new_closest_cluster != which_cluster[j])
            {
                if (which_cluster[j] != -1)
                {
                    cluster_size[which_cluster[j]] -= 1;
                    multiply_scalar((*vectors)[j], &minus_one, negative_vec);
                    add_vectors(cluster_sums[which_cluster[j]], negative_vec, cluster_sums[which_cluster[j]]);
                }
                cluster_size[new_closest_cluster] += 1;
                add_vectors(cluster_sums[new_closest_cluster], (*vectors)[j], cluster_sums[new_closest_cluster]);
                which_cluster[j] = new_closest_cluster;
                differ = 1;
            }
        }
        if (differ == 0)
        {
            break;
        }
        for (l = 0; l < k; l++)
        {
            if (cluster_size[l] != 0)
            {
                double inverse_size;
                inverse_size = (double)1 / cluster_size[l];
                multiply_scalar(cluster_sums[l], &inverse_size, (*centers)[l]);
            }
        }
    }
}