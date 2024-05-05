library(igraph)

#construct graph for adapter pairs.

adapters = c(paste("adapter",LETTERS[1:16],sep  = ""),paste("adapter",LETTERS[1:16],".rc",sep  = ""))
links = data.frame(from = adapters[1:15],to = adapters[2:16])
links2 = data.frame(from = adapters[18:32],to = adapters[17:31])
links3 = data.frame(from = "start",to = adapters[1:32])
links4 = data.frame(from = adapters[1:32],to = "end")
links = rbind(links,links2,links3,links4)
adapter_graph = graph_from_data_frame(links,directed = T)
