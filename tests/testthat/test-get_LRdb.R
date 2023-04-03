test_that("get_LRdb() gets 'human' LRdb table", {
  expect_equal(is.data.frame(get_LRdb()), TRUE)
})
