library(sequenza)

test <- sequenza.extract("288_005_small.seqz.gz")
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = "288_005", out.dir = "288_005")

test <- sequenza.extract("288_006_small.seqz.gz")
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = "288_006", out.dir = "288_006")

test <- sequenza.extract("288_008_small.seqz.gz")
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = "288_008", out.dir = "288_008")

test <- sequenza.extract("288_016_small.seqz.gz")
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = "288_016", out.dir = "288_016")