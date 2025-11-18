test_that("基础参数校验与错误处理", {
  # 创建测试数据
  test_data <- data.frame(value  = rnorm(100), category = rep(c("A","B"), 50))

  # 必需参数缺失测试
  expect_error(univar_density_barplot_plotting(plotdata = test_data), "specify the colname x")
  expect_error(univar_density_barplot_plotting(plotdata = test_data, x = "value"),
               "specify the colname group")
})
