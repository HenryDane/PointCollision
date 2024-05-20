#ifndef _PAIR_HPP
#define _PAIR_HPP

struct pair_t {
    size_t a, b;

    bool operator==(const pair_t& other) const {
        return ((a == other.a) && (b == other.b)) ||
        ((a == other.b) && (b == other.a));
    }

    bool operator<(const pair_t& other) const {
        if (a < other.a) {
            return true;
        } else if (a > other.a) {
            return false;
        } else {
            return (b < other.b);
        }
    }
};

#endif
