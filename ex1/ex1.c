#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define BUFFERSIZE  100
#define BASE_ARR_SIZE 20
#define ARR_SIZE_MULTIPLY 2

int main(int argc, char const *argv[])
{  
    int lines_count = 0;
    int arrlen = BASE_ARR_SIZE;
    char buffer[BUFFERSIZE];
    char ** lines = (char **) (malloc(arrlen * sizeof(char *)));
    assert(lines != NULL);
    while (fgets(buffer, BUFFERSIZE ,stdin) != NULL)
    {
        if (lines_count == arrlen)
        {
            arrlen *= ARR_SIZE_MULTIPLY;
            lines = (char **) (realloc(lines,arrlen * sizeof(char *)));
            assert(lines != NULL);
        }
        
        lines[lines_count] = (char *)(malloc((strlen(buffer) + 1)*sizeof(char)));
        assert(lines[lines_count] != NULL);
        strcpy(lines[lines_count],buffer);
        lines_count++;
    }
    char * data [lines_count];
    for (int i = 0; i < lines_count; i++)
    {
        data[i] = lines[i];
    }
    printf("  ");
     

    

   
    return 0;
}





