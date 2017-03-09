context("models")

test_that("photor", {
  expect_equal(ncol(photor(350)), 2)
  expect_equal(nrow(photor(350, lambda = seq(300, 700, 1))), 401)
})

test_that("photor", {
  expect_equal(ncol(photor(350, beta.band=TRUE)), 2)
  expect_equal(nrow(photor(350, lambda = seq(300, 700, 1), beta.band=TRUE)), 401)
})


test_that("logistic", {
  expect_equal(ncol(logistic(x0=350,L=60,k=0.08)), 2)
  expect_equal(nrow(logistic(x0=350,L=60,k=0.08)), 401)
})

data("bee")
data("D65")
data("Rb")

midpoint<-seq(from = 300, to = 700, 50)
W<-seq(300, 700, 1)
R<-data.frame(W)
for (i in 1:length(midpoint)) {
  R[,i+1]<-logistic(x = seq(300, 700, 1), x0=midpoint[[i]], L = 50, k=0.04)[,2]
}
names(R)[2:ncol(R)]<-midpoint

#1.CTTKmodel dichromatic
test_that("CTTKmodel", {
  expect_equal(ncol(CTTKmodel(photo=2,
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=bee)), 6)
})


#1.CTTKmodeltrichromatic
test_that("CTTKmodel", {
  expect_equal(ncol(CTTKmodel(photo=c("tri"),
                                R=R,
                                I=D65,
                                Rb=Rb,
                                C=bee)), 9)
})

test_that("CTTKmodel", {
  expect_equal(nrow(CTTKmodel(photo=c("tri"),
                                R=R,
                                I=D65,
                                Rb=Rb,
                                C=bee)), 9)
})

#2.CTTKmodel tetrachromatic
test_that("CTTKmodel", {
  expect_equal(ncol(CTTKmodel(photo=c("tetra"),
                                R=R,
                                I=D65,
                                Rb=Rb,
                                C=photor(c(350,420,490,560)))), 12)
})

test_that("CTTKmodel", {
  expect_equal(nrow(CTTKmodel(photo=c("tetra"),
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=photor(c(350,420,490,560)))), 9)
})

test_that("CTTKmodel", {
  expect_equal(ncol(CTTKmodel(photo=6,
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=photor(c(350,380,420,490,520,560)))), 18)
})

#3. EMmodel trichromatic
test_that("EMmodel", {
  expect_equal(ncol(EMmodel(photo=2,
                            R=R,
                            I=D65,
                            Rb=Rb,
                            C=bee)), 6)
})


test_that("EMmodel", {
  expect_equal(ncol(EMmodel(photo=c("tri"),
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=bee)), 9)
})

test_that("EMmodel", {
  expect_equal(nrow(EMmodel(photo=c("tri"),
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=bee)), 9)
})

#4.EMmodel tetrachromatic
test_that("EMmodel", {
  expect_equal(ncol(EMmodel(photo=c("tetra"),
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=photor(c(350,420,490,560)))), 12)
})

test_that("EMmodel", {
  expect_equal(nrow(EMmodel(photo=c("tetra"),
                              R=R,
                              I=D65,
                              Rb=Rb,
                              C=photor(c(350,420,490,560)))), 9)
})


test_that("EMmodel", {
  expect_equal(ncol(EMmodel(photo=6,
                            R=R,
                            I=D65,
                            Rb=Rb,
                            C=photor(c(350,380,420,490,520,560)))), 18)
})


#5. RNLmodel, dichromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("di"),
                            model="linear",
                            R1=R,
                            I=D65,
                            Rb=Rb,
                            C=photor(c(420,560)),
                            noise=TRUE,
                            e=c(0.1,0.05))), 13)
  expect_equal(nrow(RNLmodel(photo=c("di"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=TRUE,
                             e=c(0.1,0.05))), 9)
})

#6. RNLmodel, dichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("di"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))), 13)
  expect_equal(nrow(RNLmodel(photo=c("di"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))), 9)
})


#7. RNLmodel, dichromatic, log=T, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("di"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=TRUE,
                             e=c(0.1,0.05))), 13)
  expect_equal(nrow(RNLmodel(photo=c("di"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=TRUE,
                             e=c(0.1,0.05))), 9)
})

#8. RNLmodel, dichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("di"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))), 13)
  expect_equal(nrow(RNLmodel(photo=c("di"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))), 9)
})


#9. RNLmodel, trichromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))), 20)
  expect_equal(nrow(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))), 9)
})

#10. RNLmodel, trichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2))), 20)
  expect_equal(nrow(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1,5,2))), 9)
})


#11. RNLmodel, trichromatic, log=T, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))), 20)
  expect_equal(nrow(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))), 9)
})

#12. RNLmodel, trichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2))), 20)
  expect_equal(nrow(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2))), 9)
})


#13. RNLmodel, tetrachromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))), 27)
  expect_equal(nrow(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))), 9)
})

#14. RNLmodel, tetrachromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2))), 27)
  expect_equal(nrow(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1,5,2,2))), 9)
})


#15. RNLmodel, tetrachromatic, log=T, noise=T
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))), 27)
  expect_equal(nrow(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))), 9)
})

#16. RNLmodel, tetrachromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2))), 27)
  expect_equal(nrow(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2))), 9)
})

test_that("RNLmodel", {
  expect_equal(ncol(RNLmodel(photo=5,
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560,600)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2,2))), 26)
  expect_equal(nrow(RNLmodel(photo=5,
                             model="log",
                             R1=R,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560,600)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2,2))), 9)
})


test_that("RNLthres", {
  expect_equal(ncol(RNLthres(photo=5,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560,600)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2,2))), 3)
  expect_equal(nrow(RNLthres(photo=5,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560,600)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2,2))), 401)
})


#deltaS=0
#deltaS=0
#deltaS=0
#1.CTTKmodel trichromatic
test_that("CTTKmodel", {
  expect_equal(CTTKmodel(photo=c("tri"),
                              R=Rb,
                              I=D65,
                              Rb=Rb,
                              C=bee)$deltaS, 0)
})

#2.CTTKmodel tetrachromatic
test_that("CTTKmodel", {
  expect_equal(CTTKmodel(photo=c("tetra"),
                              R=Rb,
                              I=D65,
                              Rb=Rb,
                              C=photor(c(350,420,490,560)))$deltaS, 0)
})

#3. EMmodel trichromatic
EMRb<-data.frame(300:700, rep(7,401))
test_that("EMmodel", {
  expect_equal(round(EMmodel(photo=c("tri"),
                            R=cbind(EMRb[,1],EMRb[,2]+10),
                            I=D65,
                            Rb=EMRb,
                            C=bee)$deltaS, 2), 0)
})


#4.EMmodel tetrachromatic
test_that("EMmodel", {
  expect_equal(round(EMmodel(photo=c("tetra"),
                       R=cbind(EMRb[,1],EMRb[,2]+10),
                       I=D65,
                       Rb=EMRb,
                       C=photor(c(350,420,490,560)))$deltaS, 2), 0)
})


#5. RNLmodel, dichromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("di"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=TRUE,
                             e=c(0.1,0.05))$deltaS, 0)
})

#6. RNLmodel, dichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("di"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))$deltaS, 0)
})


#7. RNLmodel, dichromatic, log=T, noise=T
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("di"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=TRUE,
                             e=c(0.1,0.05))$deltaS, 0)
})

#8. RNLmodel, dichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("di"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(420,560)),
                             noise=F,
                             v=0.1, n=c(1,2))$deltaS, 0)
})


#9. RNLmodel, trichromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))$deltaS, 0)
})

#10. RNLmodel, trichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tri"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2))$deltaS, 0)
})


#11. RNLmodel, trichromatic, log=T, noise=T
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05))$deltaS, 0)
})

#12. RNLmodel, trichromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tri"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2))$deltaS, 0)
})


#13. RNLmodel, tetrachromatic, log=F, noise=T
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))$deltaS, 0)
})

#14. RNLmodel, tetrachromatic, log=F, noise=F
test_that("RNLmodel", {
  expect_equal(RNLmodel(photo=c("tetra"),
                             model="linear",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2))$deltaS, 0)
})


#15. RNLmodel, tetrachromatic, log=T, noise=T
test_that("RNLmodel, tetrachromatic, log=T, noise=T", {
  expect_equal(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=TRUE,
                             e=c(0.1,0.07,0.05,0.04))$deltaS, 0)
})

#16. RNLmodel, tetrachromatic, log=F, noise=F
test_that("RNLmodel, tetrachromatic, log=F, noise=F", {
  expect_equal(RNLmodel(photo=c("tetra"),
                             model="log",
                             R1=Rb,
                             I=D65,
                             Rb=Rb,
                             C=photor(c(350,420,490,560)),
                             noise=F,
                             v=0.1, n=c(1,1.5,2,2))$deltaS, 0)
})



###compare with AVICOL
r500<-logistic(x0=500,L=50,k=0.04)
test_that("CTTKmodel, trichromatic, AVICOL", {
  model<-CTTKmodel(photo="tri",
                  R=r500,
                  I=D65,
                  Rb=Rb,
                  C=bee)
  expect_equal(round(model$E1, 3), 0.301)
  expect_equal(round(model$E2, 3), 0.597)
  expect_equal(round(model$E3, 3), 0.786)
})

test_that("CTTKmodel, tetrachromatic, AVICOL", {
  model<-CTTKmodel(photo="tetra",
                   R=r500,
                   I=D65,
                   Rb=Rb,
                   C=photor(c(350,420,490,560)))
  expect_equal(round(model$E1, 3), 0.112)
  expect_equal(round(model$E2, 3), 0.520)
  expect_equal(round(model$E3, 3), 0.755)
  expect_equal(round(model$E4, 3), 0.795)
})

test_that("RNLmodel, dichromatic, AVICOL", {
  model<-RNLmodel(photo="di",
                model="log",
                R1=r500,
                I=D65,
                Rb=Rb,
                C=photor(c(420,560)),
                noise=F,
                v=0.05,
                n=c(1,1.77))
  expect_equal(round(model$e1, 3), 0.050)
  expect_equal(round(model$e2, 3), 0.038)
  expect_equal(round(model$E1_R1, 3), 0.034)
  expect_equal(round(model$E2_R1, 3), 0.588)
  expect_equal(round(model$deltaS, 3), 8.853)
})

test_that("RNLmodel, trichromatic, AVICOL", {
  model<-RNLmodel(photo="tri",
                  model="log",
                  R1=r500,
                  I=D65,
                  Rb=Rb,
                  C=bee,
                  noise=F,
                  v=0.05,
                  n=c(1,5,10))
  expect_equal(round(model$E1_R1, 3), -0.366)
  expect_equal(round(model$E2_R1, 3), 0.171)
  expect_equal(round(model$E3_R1, 3), 0.564)
  expect_equal(round(model$e1, 3), 0.050)
  expect_equal(round(model$e2, 3), 0.022)
  expect_equal(round(model$e3, 3), 0.016)
  expect_equal(round(model$deltaS, 3), 21.094)
})

test_that("RNLmodel, tetrachromatic, AVICOL", {
    model<-RNLmodel(photo="tetra",
                   model="log",
                   R1=r500,
                   I=D65,
                   Rb=Rb,
                   C=photor(c(350,420,490,560)),
                   noise=F,
                   v=0.05,
                   n=c(1, 1.9, 2.2, 2.1))
    expect_equal(round(model$E1_R1, 3), -0.900)
    expect_equal(round(model$E2_R1, 3), 0.034)
    expect_equal(round(model$E3_R1, 3), 0.488)
    expect_equal(round(model$E4_R1, 3), 0.588)
    expect_equal(round(model$e1, 3), 0.050)
    expect_equal(round(model$e2, 3), 0.036)
    expect_equal(round(model$e3, 3), 0.034)
    expect_equal(round(model$e4, 3), 0.035)
    expect_equal(round(model$deltaS, 3), 26.530)
})


test_that("RNLmodel, tetrachromatic, vismodel and tcs", {
  model<-EMmodel(photo="tetra",
                 R=r500,
                 I=D65,
                 Rb=Rb,
                 C=photor(c(350,420,490,560)))
  expect_equal(round(model$E1, 3), -4.284)
  expect_equal(round(model$E2, 3), 0.163)
  expect_equal(round(model$E3, 3), 2.322)
  expect_equal(round(model$E4, 3), 2.799)
  expect_equal(round(model$X1, 3), 1.615)
  expect_equal(round(model$X2, 3), 0.595)
  expect_equal(round(model$X3, 3), -4.534)
  expect_equal(round(model$deltaS, 3), 4.850)
})
