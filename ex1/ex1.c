#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define BUFFERSIZE 1000
#define BASE_ARR_SIZE 20
#define ARR_SIZE_MULTIPLY 2
#define COMMA ','
#define MALLOC_FOR_STRING(str) (char *)(malloc((strlen(str) + 1) * sizeof(char)))

int d = 0;
int k = 0;
int max_iter = 200;
void split_line_to_strings_nums(char *str, char ***);
void count_commas_in_str(char *str, int *d);
void turn_lines_to_vectors(int *d, int *lines_count, char ***lines, float ***vers);
int is_positive_int(const char *str);
void check_command_line_args(int *argc, char const *argv[]);

int main(int argc, char const *argv[])
{
    check_command_line_args(&argc, argv);

    int lines_count = 0;
    int arrlen = BASE_ARR_SIZE;
    char buffer[BUFFERSIZE];
    char **lines = (char **)(malloc(arrlen * sizeof(char *)));
    assert(lines != NULL);
    while (fgets(buffer, BUFFERSIZE, stdin) != NULL)
    {
        if (lines_count == arrlen)
        {
            arrlen *= ARR_SIZE_MULTIPLY;
            lines = (char **)(realloc(lines, arrlen * sizeof(char *)));
            assert(lines != NULL);
        }
        if (d == 0)
        {
            count_commas_in_str(buffer, &d);
        }

        lines[lines_count] = (char *)(malloc((strlen(buffer) + 1) * sizeof(char)));
        assert(lines[lines_count] != NULL);
        strcpy(lines[lines_count], buffer);
        lines_count++;
    }

    float **vectors;
    turn_lines_to_vectors(&d, &lines_count, &lines, &vectors);
    free(lines);
    for (int i = 0; i < lines_count; i++)
    {
        for (int k = 0; k < d; k++)
        {
            printf("%.4f ", vectors[i][k]);
        }
        printf("\n");
    }
    return 0;
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
