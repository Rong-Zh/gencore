#ifndef UMI_CLUSTER_H
#define UMI_CLUSTER_H

#include <map>
#include <cstddef>
#include <string>
#include <vector>

// Counts are original template counts within ONE existing coordinate cluster.
// Families and their members are deterministic; the first member is the root.
std::vector<std::vector<std::string>> groupUmis(
    const std::map<std::string, std::size_t>& counts, int threshold);

#endif
