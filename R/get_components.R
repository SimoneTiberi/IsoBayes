get_components = function(EC_numeric) {
    graph = lapply(seq_len(length(EC_numeric)), function(i) {
        x = EC_numeric[[i]]
        out = rbind()
        for (ii in seq_len(length(x))) {
            out = rbind(out, c(i, x[ii]))
        }
        out
    })
    graph = do.call("rbind", graph)

    components = list()
    while (nrow(graph) > 0) {
        nodes = matrix(graph[1, ], ncol = 2)
        graph = matrix(graph[-1, ], ncol = 2)
        new_nodes = c( which(graph[, 1] %in% unique(nodes[, 1])), which(graph[, 2] %in% unique(nodes[, 2])))
        check = length(new_nodes) > 0
        while (check) {
            nodes = rbind(nodes, graph[new_nodes, ])
            graph = matrix(graph[-new_nodes, ], ncol = 2)
            new_nodes = c(which(graph[, 1] %in% unique(nodes[, 1])), which(graph[, 2] %in% unique(nodes[, 2])))
            check = length(new_nodes) > 0
        }
        components = append(components,
                            list(list(peptides = unique(nodes[, 1]),
                                      proteins = unique(nodes[, 2]))))
    }
    components
}