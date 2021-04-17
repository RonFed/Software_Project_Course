#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

int a,b,letters_needed_a;
int num_base_ten = 0;
void base_read(char t),read_num();
int legal_char(int c), char_val(int c);
int convert_to_base_b();
char char_from_b(int n);


int main()
{
    printf("Please enter the numbers base:\n");
    base_read('i');
    letters_needed_a = (a > 10);
    printf("Please enter the desired base:\n");
    base_read('d');
    printf("Please enter a number in base %d:\n",a);
    getchar();
    read_num();
    printf("The result is: ");
    convert_to_base_b();
    printf("\n");
    return 0;
}

void read_num()
{
    int d = 0;
    char c;
    while ((c = getchar()) != EOF && c != '\n') {
        if (!legal_char(c)) {
            printf("Invalid number!\n");
            exit(0);
        }
        d = char_val(c);
        num_base_ten *= a;
        num_base_ten += d;     
    }
    if (num_base_ten == 0) {
        printf("The result is: 0");
        exit(0);
    }
}


int convert_to_base_b() {
    int d;
    if (num_base_ten == 0) {
        return 0;
    }
    d = num_base_ten%b;
    num_base_ten /= b;
    convert_to_base_b();
    printf("%c",char_from_b(d));
    return 0;

}

char char_from_b(int n) {
    if (n <= 9) {
        return '0' + n;
    } else {
        return 'A' + n - 10;
    }
}


int char_val(int c)
{
    if (c <= '9') {
        return c - '0';
    }
    return c - 'A' + 10; 
}

int legal_char(int c)
{
    if (letters_needed_a) {
        if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'A' + a -11))) {
            return 0;
        }
    } else {
        if (c < '0' || c >=  '0' + a) {
            return 0;
        }
    }
    return 1;
}

void base_read(char t) 
{
    if (t == 'i') {
        scanf("%d", &a);
        if (a < 2 || a > 16) {
            printf("Invalid input base\n");
            exit(0);
        }
    }
    if (t == 'd') {
        scanf("%d", &b);
        if (b < 2 || b > 16) {
            printf("Invalid desired base\n");
            exit(0);
        }
    }
    return;
}
