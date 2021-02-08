#pragma once

#include <fstream>
#include <iostream>
#include <chrono>

class FileLog {
public:
    static FileLog& getInstance(const std::string &logFile = "") {
        static FileLog instance(logFile);
        return instance;
    }

    static void setUpFileLog(const std::string &logFile) {
        getInstance(logFile);
    }
    
    static void setPrefix(const std::string &prefix) {
        getInstance().prefix = prefix;
    }
    
    std::string prefix;
    std::ofstream stream;
    bool isInited = false;
    
private:
    FileLog(const std::string &logFile) : stream(logFile) {
        if (logFile != "") {
            isInited = true;
        }
    }

public:
    FileLog(FileLog const&) = delete;
    void operator=(FileLog const&) = delete;
};

template <typename T>
std::ofstream& operator<<(FileLog& log, const T &t) {
    log.stream << log.prefix << t;
    if (!log.stream) {
        std::cout << "flog: " << t << std::endl;
    }
    return log.stream;
}

static inline FileLog& flog() {
    return FileLog::getInstance();
}

static inline double secondsBetween(clock_t fromTime, clock_t toTime) {
    return ((double)(toTime - fromTime) / CLOCKS_PER_SEC);
}
    
static inline void logTime(const std::string &description = "", clock_t fromTime = clock_t()) {
    static clock_t lastTime = clock();
    clock_t currentTime = clock();
    if (description.empty()) {
        auto seconds = fromTime == clock_t() ? secondsBetween(lastTime, currentTime) : secondsBetween(fromTime, currentTime);
        if (flog().isInited) {
            flog() << description << " time: " << seconds << "s\n";
        } else {
            std::cout << description << " time: " << seconds << "s\n";
        }
    }
    lastTime = currentTime;
}

struct LogGuard {
    using Clock = std::chrono::steady_clock;

    const std::string id_;
    const Clock::time_point start_time_ = Clock::now();

    LogGuard(const std::string &id) : id_(id) {}

    ~LogGuard() {
        const auto end_time = Clock::now();
        const auto dur = end_time - start_time_;
        auto mills = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        std::cout << id_ << " time " << mills / 1000. << "s\n";
    }
};
