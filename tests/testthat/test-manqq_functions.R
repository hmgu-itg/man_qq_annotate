

test_that("refine_data", {
  ## Case 1
  # p-value column read as.character due to
  # a non-numeric value (e.g. 'nan') needs to be
  # converted to numeric and filtered.
  data = data.table::data.table(
      p = c('0.1', '0.1', 'nan')
  )
  refined_data = refine_data(data)
  
  expected = data.table::data.table(
      p = c(0.1, 0.1)
  )
  
  expect_true(all.equal(refined_data, expected))
})