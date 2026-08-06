test_that("df_to_graph returns an igraph object", {
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"), stringsAsFactors = FALSE)
  g <- df_to_graph(edges, gc = FALSE)
  expect_true(igraph::is_igraph(g))
})

test_that("df_to_graph without nodelist has correct vertex count", {
  edges <- data.frame(from = c("A", "B"), to = c("B", "C"), stringsAsFactors = FALSE)
  g <- df_to_graph(edges, gc = FALSE)
  expect_equal(igraph::vcount(g), 3L)
  expect_equal(igraph::ecount(g), 2L)
})

test_that("df_to_graph with nodelist filters correctly", {
  edges <- data.frame(
    from = c("A", "B", "D"),
    to   = c("B", "C", "E"),
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(name = c("A", "B", "C"), stringsAsFactors = FALSE)
  g <- df_to_graph(edges, nodelist = nodes, gc = FALSE)
  # Edges A->B and B->C remain (both endpoints in nodelist); D->E is dropped
  expect_equal(igraph::vcount(g), 3L)
  expect_equal(igraph::ecount(g), 2L)
})

test_that("df_to_graph gc=TRUE returns the giant component", {
  # Two disconnected components: A-B-C and D-E
  edges <- data.frame(
    from = c("A", "B", "D"),
    to   = c("B", "C", "E"),
    stringsAsFactors = FALSE
  )
  # Suppress igraph deprecation warning from package internals
  g <- suppressWarnings(df_to_graph(edges, gc = TRUE))
  # Giant component is A-B-C (3 nodes)
  expect_equal(igraph::vcount(g), 3L)
})

test_that("df_to_graph gc=FALSE returns all nodes when graph is connected", {
  edges <- data.frame(from = c("A", "B", "C"), to = c("B", "C", "D"), stringsAsFactors = FALSE)
  g_gc  <- suppressWarnings(df_to_graph(edges, gc = TRUE))
  g_all <- df_to_graph(edges, gc = FALSE)
  expect_equal(igraph::vcount(g_gc), igraph::vcount(g_all))
})

test_that("df_to_graph with null nodelist and gc=FALSE preserves all edges", {
  from <- c("1", "2", "4")
  to   <- c("2", "3", "5")
  edges <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
  g <- df_to_graph(edges, gc = FALSE)
  expect_equal(igraph::ecount(g), 3L)
})

test_that("df_to_graph example from documentation works", {
  from  <- c("1", "2", "4")
  to    <- c("2", "3", "5")
  edges <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
  nodes <- data.frame(name = c("1", "2", "3", "4", "5"), stringsAsFactors = FALSE)
  g <- suppressWarnings(df_to_graph(edges, nodes))
  expect_true(igraph::is_igraph(g))
})
