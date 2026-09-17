#include "umicluster.h"

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

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

std::vector<std::vector<std::string>> groupUmis(
    const std::map<std::string, size_t>& counts, int threshold) {
    if(threshold < 0)
        throw std::invalid_argument("UMI distance threshold cannot be negative");
    struct Node { std::string umi; size_t count; };
    std::vector<Node> nodes;
    nodes.reserve(counts.size());
    for(const auto& [umi, count] : counts) {
        if(count > 0)
            nodes.push_back({umi, count});
    }
    std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a.count != b.count ? a.count > b.count : a.umi < b.umi;
    });

    std::vector<std::vector<std::string>> families;
    families.reserve(nodes.size());
    // Small clusters are cheaper to scan. Larger distance-one clusters can
    // find all neighbors by substituting one base and probing an index.
    bool indexed = threshold == 1 && nodes.size() > 512;
    for(const auto& node : nodes) {
        if(node.umi.find_first_not_of("ACGTN_") != std::string::npos)
            indexed = false;
    }
    std::unordered_map<std::string_view, size_t> index;
    if(indexed) {
        index.reserve(nodes.size());
        for(size_t i = 0; i < nodes.size(); ++i)
            index.emplace(nodes[i].umi, i);
    }
    std::vector<size_t> neighbors;
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
                const auto assign = [&](size_t child) {
                    if(assigned[child])
                        return;
                    // Equivalent to parent >= 2*child-1, without overflow.
                    if(nodes[child].count > parent.count / 2 + parent.count % 2)
                        return;
                    assigned[child] = true;
                    family.push_back(nodes[child].umi);
                    pending.push_back(child);
                };
                if(indexed) {
                    neighbors.clear();
                    std::string candidate = parent.umi;
                    for(size_t position = 0; position < candidate.size(); ++position) {
                        const char original = candidate[position];
                        if(original == '_')
                            continue;
                        for(char base : std::string_view("ACGTN")) {
                            if(base == original)
                                continue;
                            candidate[position] = base;
                            const auto found = index.find(candidate);
                            if(found != index.end() && !assigned[found->second])
                                neighbors.push_back(found->second);
                        }
                        candidate[position] = original;
                    }
                    // Preserve the old count/lexical traversal order exactly.
                    std::sort(neighbors.begin(), neighbors.end());
                    for(size_t child : neighbors)
                        assign(child);
                } else {
                    for(size_t child = 0; child < nodes.size(); ++child) {
                        if(!assigned[child] &&
                           nodes[child].count <= parent.count / 2 + parent.count % 2 &&
                           withinDistance(parent.umi, nodes[child].umi, threshold))
                            assign(child);
                    }
                }
            }
        }
        families.push_back(std::move(family));
    }
    return families;
}
