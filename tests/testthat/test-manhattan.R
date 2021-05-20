# Tests for functions in the manhattan.R script

test_that("process_overlap_restr", {
  ## Case 1
  # One protein coding gene overlap
  one_external_name = data.frame(
                      biotype = c('protein_coding'),
                      external_name = c('AAA')
                      )
  gene_name = process_overlap_restr(one_external_name)
  expect_equal(gene_name, 'AAA')


  ## Case 2
  # Multiple protein coding gene overlap
  external_names = data.frame(
                        biotype = c('protein_coding', 'protein_coding'),
                        external_name = c('AAA', 'BBB')
                        )

  gene_name = process_overlap_restr(external_names)
  expect_equal(gene_name, 'AAA,BBB')


  ## Case 3
  # A protein coding, novel gene overlap but gene has no external name.
  gene_id = 'ENSG00000254469'
  no_external_name = data.frame(
                      biotype = c('protein_coding'),
                      id = c(gene_id)
                      )
  # Output Ensembl gene stable ID instead because this gene has no external_name
  gene_name = process_overlap_restr(no_external_name)
  expect_equal(gene_name, gene_id)


  ## Case 4
  # A protein coding, novel gene overlap but gene has no external name.
  no_external_name = data.frame(
                      biotype = c('protein_coding', 'protein_coding'),
                      id = c('ENSG00000254469', 'ENSG00000300000')
                      )
  # Output Ensembl gene stable ID instead because this gene has no external_name
  gene_name = process_overlap_restr(no_external_name)
  expect_equal(gene_name, 'ENSG00000254469,ENSG00000300000')

})


test_that("query_ensembl_gene_overlap", {
  ## Case 1
  # There's no gene at position chr1:1, so we expect an empty data.frame
  restr = query_ensembl_gene_overlap(1, 1, 1)

  expect_s3_class(restr, "data.frame")
  expect_equal(nrow(restr), 0)

  ## Case 2
  # The position overlaps the famous NODAL gene.
  restr = query_ensembl_gene_overlap(10, 70431936, 70431937)
  expect_s3_class(restr, "data.frame")
  expect_equal(nrow(restr), 1)
})