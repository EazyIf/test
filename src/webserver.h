#ifndef WEBSERVER_H
#define WEBSERVER_H

#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstring>
#include <cstdio>
#include <string>
#include <functional>
#include <map>

struct HttpRequest {
    std::string method;
    std::string path;
    std::string body;
};

class WebServer {
public:
    using Handler = std::function<std::string(const HttpRequest&)>;

    bool start(int port) {
        listen_fd_ = socket(AF_INET, SOCK_STREAM, 0);
        if (listen_fd_ < 0) return false;

        int opt = 1;
        setsockopt(listen_fd_, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY;
        addr.sin_port = htons(port);

        if (bind(listen_fd_, (sockaddr*)&addr, sizeof(addr)) < 0) {
            close(listen_fd_); return false;
        }
        if (listen(listen_fd_, 8) < 0) {
            close(listen_fd_); return false;
        }

        // Non-blocking
        int flags = fcntl(listen_fd_, F_GETFL, 0);
        fcntl(listen_fd_, F_SETFL, flags | O_NONBLOCK);

        port_ = port;
        return true;
    }

    void route(const std::string& path, Handler handler) {
        routes_[path] = handler;
    }

    // Call once per frame — handles at most one pending connection
    void poll() {
        sockaddr_in client_addr{};
        socklen_t client_len = sizeof(client_addr);
        int client_fd = accept(listen_fd_, (sockaddr*)&client_addr, &client_len);
        if (client_fd < 0) return;

        // Set a short timeout for reading
        struct timeval tv{0, 100000}; // 100ms
        setsockopt(client_fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));

        char buf[4096];
        int n = recv(client_fd, buf, sizeof(buf) - 1, 0);
        if (n <= 0) { close(client_fd); return; }
        buf[n] = '\0';

        HttpRequest req = parse_request(buf, n);

        std::string response;
        auto it = routes_.find(req.path);
        if (it != routes_.end()) {
            response = it->second(req);
        } else {
            response = "HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n";
        }

        // Send response (may need multiple writes for large responses)
        const char* data = response.data();
        size_t remaining = response.size();
        while (remaining > 0) {
            ssize_t sent = send(client_fd, data, remaining, MSG_NOSIGNAL);
            if (sent <= 0) break;
            data += sent;
            remaining -= sent;
        }

        close(client_fd);
    }

    void stop() {
        if (listen_fd_ >= 0) { close(listen_fd_); listen_fd_ = -1; }
    }

    int port() const { return port_; }

    ~WebServer() { stop(); }

private:
    int listen_fd_ = -1;
    int port_ = 0;
    std::map<std::string, Handler> routes_;

    HttpRequest parse_request(const char* buf, int len) {
        HttpRequest req;
        std::string raw(buf, len);

        // Parse method
        size_t sp1 = raw.find(' ');
        if (sp1 == std::string::npos) return req;
        req.method = raw.substr(0, sp1);

        // Parse path
        size_t sp2 = raw.find(' ', sp1 + 1);
        if (sp2 == std::string::npos) return req;
        req.path = raw.substr(sp1 + 1, sp2 - sp1 - 1);

        // Strip query string
        size_t q = req.path.find('?');
        if (q != std::string::npos) req.path = req.path.substr(0, q);

        // Parse body (after \r\n\r\n)
        size_t body_start = raw.find("\r\n\r\n");
        if (body_start != std::string::npos) {
            req.body = raw.substr(body_start + 4);
        }

        return req;
    }
};

static std::string http_response(const std::string& content_type, const std::string& body) {
    std::string resp = "HTTP/1.1 200 OK\r\n";
    resp += "Content-Type: " + content_type + "\r\n";
    resp += "Content-Length: " + std::to_string(body.size()) + "\r\n";
    resp += "Access-Control-Allow-Origin: *\r\n";
    resp += "Cache-Control: no-cache\r\n";
    resp += "\r\n";
    resp += body;
    return resp;
}

#endif
