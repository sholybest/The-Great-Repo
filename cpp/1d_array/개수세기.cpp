#include <iostream>

int main(){

    int num1, num2, num3;
    
    std::cin >> num1;
    int arr[num1];

    for(int i=0; i<num1; i++){
        std::cin >> num2;
        arr[i] = num2;
    }

    std::cin >> num3;

    int cnt = 0;
    for(int i=0; i<num1; i++){
        cnt += (int)(arr[i] == num3);
    }

    std::cout << cnt << std::endl;

    return 0;
}