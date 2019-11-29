#ifndef COLOR_H
#define COLOR_H
#include <cinttypes>

class Color {
public:
    uint8_t r, g, b;

    Color(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0);
};

#endif // COLOR_H
