## Build essential functions
# Lea-Coulson function
fun.m.rate.lc = function(m, r) {
  return(abs(1.24 - ((r/m) - log(m))))
}

# "Average rate" function
fun.m.rate.avg = function(x, N , C, r) {
  return(abs((x * N * log(x * N * C)) - r))
}


# Find solution to Lea-Coulson function
solve.m.rate.lc = function(mutants, N) {
  r = median(mutants[,1])
  opt = optimize(f=fun.m.rate.lc, r=r, lower = 0, upper = 1E15, tol=1E-15)
  return(c(rate = format(opt$minimum/(N),scientific=TRUE), opt.prec = format(opt$objective,scientific=TRUE)))
}

# Find solution to Average Rate function
solve.m.rate.avg = function(mutants, N) {
  r = mean(mutants[,1])
  C = length(mutants[,1])
  opt = optimize(f=fun.m.rate.avg, N=N, C=C, r=r, lower = 0, upper = 1E15, tol=1E-15)
  return(c(rate = format(opt$minimum,scientific=TRUE), opt.prec = format(opt$objective,scientific=TRUE)))
}

# Find mutation rate using proportion of 0 mutant cultures
m.rate.zeros = function(mutants,N) {
  p0 = sum(mutants[,1] == 0)/nrow(mutants)
  return(c(rate=format(-log(p0)/N,scientific=TRUE),opt.prec=NA))
}

## Imports
library(shiny)
library(readxl)
library(tools)

shinyServer(function(input, output) {
  output$out <- renderTable({
      
    message(input$mutants$datapath)
    
    if (!is.null(input$mutants$datapath))
    {
      if(file_ext(input$mutants$name) == 'xlsx')
      {
            if(is.null(input$mutants))
                return(NULL)

            ext = 'xlsx'
            if (file.exists(input$mutants$datapath))
            {
                file.rename(input$mutants$datapath, paste(input$mutants$datapath, ext, sep="."))
            }
            mutants = read_excel(paste(input$mutants$datapath, ext, sep="."), col_names=TRUE)['mutants']
      } else {
          mutants = read.csv(input$mutants$datapath, header=TRUE)['mutants']
      }
      N = input$N
      
      table = data.frame(Zeros = m.rate.zeros(mutants,N), LeaCoulson = solve.m.rate.lc(mutants,N), Average= solve.m.rate.avg(mutants, N),row.names=c('Mutation Rate', 'Precision'))
      print(table)
      #print(mutants)
    }
  },rownames = TRUE)
})
