#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define BUFFERSIZE 1000
#define BASE_ARR_SIZE 20
#define ARR_SIZE_MULTIPLY 2
#define COMMA ','
#define MALLOC_FOR_STRING(str) (char *)(malloc((strlen(str) + 1) * sizeof(char)))
void updateS(float m[], int s[], int x);

int d = 0;
int k = 0;
float minus_one = -1;
int max_iter = 200;
int lines_count;
void read_text_to_lines(char *** lines);
void split_line_to_strings_nums(char *str, char ***);
void count_commas_in_str(char *str, int *d);
void turn_lines_to_vectors(int *d, int *lines_count, char ***lines, float ***vers);
int is_positive_int(const char *str);
void check_command_line_args(int *argc, char const *argv[]);
void distacne(float **a, float **b, float *result);
void closest_cluster_index(float *vector, float **centers, int *index_result);
void multiply_scalar(float *vector, float *scalar, float *result);
void add_vectors(float *a, float *b, float *result);
void print_final(float **mat);
void malloc_for_mat_and_set_zeros(float ***mat, int lines, int cols);

int main(int argc, char const *argv[])
{
    check_command_line_args(&argc, argv);

    char **lines;
    read_text_to_lines(&lines);

    float **vectors;
    turn_lines_to_vectors(&d, &lines_count, &lines, &vectors);
    free(lines);

    float **centers;
    malloc_for_mat_and_set_zeros(&centers, k, d);
    float **cluster_sums;
    malloc_for_mat_and_set_zeros(&cluster_sums, k, d);
    int cluster_size[k];
    int which_cluster[lines_count];
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < d; j++)
        {
            centers[i][j] = vectors[i][j];
        }
        cluster_size[i] = 0;
    }
    for (int i = 0; i < lines_count; i++)
    {
        which_cluster[i] = -1;
    }
    for (int i = 0; i < max_iter; i++)
    {
        int differ = 0;
        for (int j = 0; j < lines_count; j++)
        {
            int new_closest_cluster;
            closest_cluster_index(vectors[j], centers, &new_closest_cluster);
            if (new_closest_cluster != which_cluster[j])
            {
                if (which_cluster[j] != -1)
                {
                    cluster_size[which_cluster[j]] -= 1;
                    float *negative_vec;
                    negative_vec = (float *)malloc(d * sizeof(float));
                    multiply_scalar(vectors[j], &minus_one, negative_vec);
                    add_vectors(cluster_sums[which_cluster[j]], negative_vec, cluster_sums[which_cluster[j]]);
                }
                cluster_size[new_closest_cluster] += 1;
                add_vectors(cluster_sums[new_closest_cluster], vectors[j], cluster_sums[new_closest_cluster]);
                which_cluster[j] = new_closest_cluster;
                differ = 1;
            }
        }
        if (differ == 0)
        {
            break;
        }
        for (int l = 0; l < k; l++)
        {
            if (cluster_size[l] != 0)
            {
                float inverse_size;
                inverse_size = (float) 1 / cluster_size[l];
                multiply_scalar(cluster_sums[l], &inverse_size, centers[l]);
            }
        }
    }
    print_final(centers);

    return 0;
}

void read_text_to_lines(char *** lines) {
    lines_count = 0;
    int arrlen = BASE_ARR_SIZE;
    char buffer[BUFFERSIZE];
    *lines = (char **)(malloc(arrlen * sizeof(char *)));
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
    for (int i = 0; i < strlen(str); i++)
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

    *result = (char **)malloc(d * sizeof(char *));
    char *num_str = strtok(str, ",");
    int j = 0;
    while (num_str)
    {
        (*result)[j] = MALLOC_FOR_STRING(num_str);
        strcpy((*result)[j], num_str);
        j++;
        num_str = strtok(NULL, ",");
    }
}

void turn_lines_to_vectors(int *d, int *lines_count, char ***lines, float ***vers)
{
    (*vers) = (float **)(malloc((*lines_count) * sizeof(float *)));
    char *data[*lines_count];
    char **current_line = NULL;
    for (int i = 0; i < *lines_count; i++)
    {
        (*vers)[i] = (float *)(malloc((*d) * sizeof(float)));
        data[i] = (*lines)[i];
        data[i][strcspn(data[i], "\n")] = '\0';
        split_line_to_strings_nums(data[i], &current_line);
        for (int j = 0; j < *d; j++)
        {
            (*vers)[i][j] = (float)atof(current_line[j]);
        }
    }
}

int is_positive_int(const char *str)
{
    char digit = str[0];
    if (digit <= '0' || digit > '9')
    {
        return 0;
    }
    for (int i = 1; i < strlen(str); i++)
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

void distacne(float **a, float **b, float *result)
{
    float sum = 0;
    for (int i = 0; i < d; i++)
    {
        sum += ((*a)[i] - (*b)[i]) * ((*a)[i] - (*b)[i]);
    }
    *result = sum;
}

void closest_cluster_index(float *vector, float **centers, int *index_result)
{
    int min_dist_index = 0;
    float min_dist;
    distacne(&vector, &centers[0], &min_dist);
    for (int i = 0; i < k; i++)
    {
        float current_dist;
        distacne(&vector, &centers[i], &current_dist);
        if (current_dist < min_dist)
        {
            min_dist = current_dist;
            min_dist_index = i;
        }
    }
    *index_result = min_dist_index;
}

void multiply_scalar(float *vector, float *scalar, float *result)
{
    for (int i = 0; i < d; i++)
    {
        result[i] = vector[i] * (*scalar);
    }
}

void add_vectors(float *a, float *b, float *result)
{
    for (int i = 0; i < d; i++)
    {
        result[i] = a[i] + b[i];
    }
}

void print_final(float **mat)
{
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < d - 1; j++)
        {
            printf("%.4f,", mat[i][j]);
        }
        printf("%.4f", mat[i][d - 1]);
        printf("\n");
    }
}

void malloc_for_mat_and_set_zeros(float ***mat, int lines, int cols)
{
    *mat = (float **)malloc(lines * sizeof(float *));
    for (int i = 0; i < lines; i++)
    {
        (*mat)[i] = (float *)malloc(cols * sizeof(float));
        for (int j = 0; j < cols; j++)
        {
            (*mat)[i][j] = 0;
        }
        
    }
}