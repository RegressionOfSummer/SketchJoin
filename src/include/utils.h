#pragma once
#include<iostream>
#include<vector>
using namespace std;

bool isPrime(int num) {
    if (num <= 1) {
        return false;
    }
    if (num <= 3) {
        return true;
    }
    if (num % 2 == 0 || num % 3 == 0) {
        return false;
    }
    for (int i = 5; i * i <= num; i += 6) {
        if (num % i == 0 || num % (i + 2) == 0) {
            return false;
        }
    }
    return true;
}

std::vector<int> closestPrimes(int num, int n) {
    std::vector<int> primes;
    int i = num;
    while (primes.size() < n) {
        if (isPrime(i)) {
            primes.push_back(i);
        }
        --i;
    }

    return primes;
}

unsigned int M_Random_Generate()
{
  	unsigned int x = rand();
  	unsigned int h = rand();

  	return x ^ ((h & 1) << 31);
}