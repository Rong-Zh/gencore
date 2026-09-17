// Frozen pre-optimization implementation from commit 8a96b90 for differential tests.
#include "umicluster.h"

#include <algorithm>
#include <stdexcept>

namespace {
bool withinDistance(const std::string& a, const std::string& b,
                    int threshold) {
    if(a == b)
        return true;
    // Missing UMIs never supply evidence for an observed UMI.
    if(a.empty() || b.empty() || a.size() != b.size())
        return false;
    size_t differences = 0;
    for(size_t i = 0; i < a.size(); ++i) {
        // The separator defines dual-UMI structure, not a mutable base.
        if((a[i] == '_') != (b[i] == '_'))
            return false;
        differences += a[i] != b[i];
        if(differences > static_cast<size_t>(threshold))
            return false;
    }
    return differences <= static_cast<size_t>(threshold);
}
}

std::vector<std::vector<std::string>> quadraticGroupUmis(
    const std::map<std::string, size_t>& counts, int threshold) {
    if(threshold < 0)
        throw std::invalid_argument("UMI distance threshold cannot be negative");
    struct Node { std::string umi; size_t count; };
    std::vector<Node> nodes;
    for(const auto& [umi, count] : counts) {
        if(count > 0)
            nodes.push_back({umi, count});
    }
    std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a.count != b.count ? a.count > b.count : a.umi < b.umi;
    });

    std::vector<std::vector<std::string>> families;
    std::vector<bool> assigned(nodes.size(), false);
    for(size_t root = 0; root < nodes.size(); ++root) {
        if(assigned[root])
            continue;
        assigned[root] = true;
        std::vector<size_t> pending{root};
        std::vector<std::string> family{nodes[root].umi};
        if(threshold > 0) {
            // Discover outgoing edges on demand: no quadratic adjacency matrix.
            for(size_t cursor = 0; cursor < pending.size(); ++cursor) {
                const auto& parent = nodes[pending[cursor]];
                for(size_t child = 0; child < nodes.size(); ++child) {
                    if(assigned[child])
                        continue;
                    // Equivalent to parent >= 2*child-1, without overflow.
                    if(nodes[child].count > parent.count / 2 + parent.count % 2)
                        continue;
                    if(!withinDistance(parent.umi, nodes[child].umi, threshold))
                        continue;
                    assigned[child] = true;
                    family.push_back(nodes[child].umi);
                    pending.push_back(child);
                }
            }
        }
        families.push_back(std::move(family));
    }
    return families;
}
