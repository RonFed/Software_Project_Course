#include <stdio.h>
#include <assert.h>
#include <math.h>

int isLegit(char c, int a);
int isDigit(char c);
int isLetterDigit(char c);


int main(void)
{
	int a;
	int b;
	int ans = 0;
	int i = 0;
	int digit = 0;
	char c;
	int k;
	int* p=&k;
	int val;
	printf("Please enter the numbers base:\n");
	scanf("%d", &a);
	if (a > 16 || a < 2) {
		printf("Invalid input base\n");
		assert(!(a > 16 || a < 2));
	}
	printf("Please enter the desired base:\n");
	scanf("%d", &b);
	if (b > 16 || b < 2) {
		printf("Invalid desired base\n");
		assert(!(b > 16 || b < 2));
	}
	printf("Please enter a number in base %d:\n", a);
	c = getchar();
	while (((c = getchar()) != EOF && c != '\n')) {
		if (isLegit(c, a) == 0) {
			printf("Invalid number!\n");
			assert(0);
		}

		if (isDigit(c) == 1) {
			digit = (int)c - 48;
		}
		if (isLetterDigit(c) == 1) {
			digit = (int)c - 97 + 10;
		}
		p++;
		int* val_p = 0;
		val_p = &digit;
		*p = *val_p;
		printf("The p is: %d:%d\n", *p,p);
		i++;
	}
	i--;
	int ii = i;
	for (; i >= 0; i--) {
		printf("The p is: %d:%d\n", *p, p);
		ans += (*p) * (pow(a, ii-i));
		p--;
	}
	printf("The result is: %d:\n", ans);
	return 0;
}

int isLegit(char c, int a)
{
	if (a < 10) {
		if ((int)c < 48 || (int)c >(47 + a)) {
			return 0;
		}
	}
	else {
		if ((int)c < 48 || ((int)c >57 && (int)c < 97) || (int)c >96+a-10){
			return 0;
		}
	}
}

int isDigit(char c)
{
	if ((int)c < 48 || (int)c > 57) {
		return 0;
	}
	return 1;
}
int isLetterDigit(char c)
{
	if ((int)c < 97 || (int)c > 102) {
		return 0;
	}
	return 1;
}