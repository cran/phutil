using("tinysnapshot")

opts <- options(cli.width = 80)

# Create test data
test_persistence <- list(
  as_persistence(list(matrix(1:4, 2, 2))),
  as_persistence(list(matrix(1:4, 2, 2)))
)

# Test valid input
expect_silent(as_persistence_set(test_persistence))
expect_true(inherits(as_persistence_set(test_persistence), "persistence_set"))

# Test error for non-list input
expect_error(as_persistence_set(1))
expect_error(as_persistence_set("string"))

# Test error for empty list
expect_error(as_persistence_set(list()))

# Test error for non-persistence elements
expect_error(as_persistence_set(list(1, 2)))
expect_error(as_persistence_set(list(test_persistence[[1]], "not_persistence")))

# Test class preservation
result <- as_persistence_set(test_persistence)
expect_true("list" %in% class(result))
expect_true("persistence_set" %in% class(result))

# Test format and print methods
expect_snapshot_print(result, label = "print-persistence-set-class")

options(opts)
