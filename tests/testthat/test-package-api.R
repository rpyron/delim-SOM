test_that("public delimSOM functions are exported", {

  expected_exports <- c(
    "train.SOM",
    "clustering.SOM",
    "plot.structure.SOM",
    "plot.learning.SOM",
    "plot.layer.distance.scale.SOM",
    "plot.K.SOM",
    "plot.model.SOM",
    "plot.map.SOM",
    "plot.variable.importance.SOM",
    "process.SNP.data.SOM",
    "plot.layer.importance.varimp.SOM",
    "plot.layer.importance.leaveoneout.SOM",
    "make.cols.binary.SOM",
    "remove.lowCV.multicollinearity.SOM"
  )

  package_exports <- getNamespaceExports("delimSOM")

  expect_true(
    all(expected_exports %in% package_exports)
  )
})



test_that("public delimSOM exports are functions", {

  expected_exports <- c(
    "train.SOM",
    "clustering.SOM",
    "plot.structure.SOM",
    "plot.learning.SOM",
    "plot.layer.distance.scale.SOM",
    "plot.K.SOM",
    "plot.model.SOM",
    "plot.map.SOM",
    "plot.variable.importance.SOM",
    "process.SNP.data.SOM",
    "plot.layer.importance.varimp.SOM",
    "plot.layer.importance.leaveoneout.SOM",
    "make.cols.binary.SOM",
    "remove.lowCV.multicollinearity.SOM"
  )

  package_namespace <- asNamespace("delimSOM")

  function_check <- vapply(
    expected_exports,
    function(fun) {
      is.function(
        get(
          fun,
          envir = package_namespace
        )
      )
    },
    FUN.VALUE = logical(1)
  )

  expect_true(
    all(function_check)
  )
})
