#include <iostream>
#include <iomanip>

int main() {
    int div = 2;
    float test = (float) 1 / (float) div;
    std::cout << std::fixed;
    std::cout << std::setprecision(4);
    std::cout << "Test is: " << test << "\n";

    return 0;
}