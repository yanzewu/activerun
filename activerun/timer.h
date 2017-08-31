#pragma once

#include <vector>
#include <time.h>

class Timer {
public:

    void init(const std::vector<const char*> names) {
        this->names = names;
        acc_time.resize(names.size());
        last_time = 0;
    }

    void set_checkpoint(int id) {
        int64_t total_time = clock();
        if(id >= 0)acc_time[id] += total_time - last_time;
        last_time = total_time;
    }
    
    void print(int check_circles, int total_circles) {
        int64_t total_time = clock();
        printf("\n\nTotal time: %.2fs\n", (double)total_time / CLOCKS_PER_SEC);
        printf("Using %d circles for estimation:\n", check_circles);
        double scale = (double)total_circles / check_circles;
        for (size_t i = 0; i < names.size(); i++) {
            printf("%s time:\t\t\t %.2f%%\n", names[i], 100.0 * scale * acc_time[i] / total_time);
        }
    }

private:

    int64_t last_time;

    std::vector<int64_t> acc_time;
    std::vector<const char*> names;
};