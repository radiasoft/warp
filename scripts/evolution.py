evolution_version = "$Id: evolution.py,v 1.1 2000/10/16 18:34:19 dave Exp $"
###########################################################################
# Performs global optimization using genetic evolution algorithms.
# Algorithm taken from Dr Dobb's Journal, April 1997, K Price, R. Storn
###########################################################################

<< "Differential Evolution"
<< "Input:"
<< "  npop = size of population (must be greater than 3)"
<< "  nparams = number of parameters"
<< "  crossover = fraction of crossovers, in range [0,1)"
<< "  f = differential factor, in range (0,1.2]"
<< "  evaluate(params) is function, given a set of parameters, returns a score"
<< "Output:"
<< "  best_params = array hold parameters which give the lowest score"
<< "The function evolve_init can be called before evolve to initialize"
<< "the population or it can be done by hand."

if(~exists("crossover")) chameleon crossover = .5
if (crossover < 0. | crossover > 1.) \
  << "Warning: crossover outside of the range [0,1)"

if(~exists("f")) chameleon f = .7
if (f < 0. | f > 1.2) \
  << "Warning: differential factor f outside of the range (0,1.2]"

if(~exists("npop")) chameleon npop = 5
if (npop < 4) then
  << "Error: number of populations, npop, must be greater than 3"
  kaboom(0)
endif

if(~exists("nparams")) chameleon nparams = 2

double trial(nparams)
double x1(npop,nparams), x2(npop,nparams)
double cost(npop)

# --- Function to return best set of parameters so far
function best_params
  return x1(mnx(cost),)
endf

# --- Function to initialize the population
# --- Default is to pick parameters randomly distributed by 1% about a
# --- base sample set of parameters.
function evolve_init(sample;deltas)
  default(deltas) = ones(shape(sample))/100.  # deltas is same size as sample
  integer i
  x1(1,) = sample
  cost(1) = evaluate(sample)
  do i=2,npop
    x1(i,) = sample*(1.+2.*(ranf(sample)-.5)*deltas)
    cost(i) = evaluate(x1(i,))
  enddo
endf

# --- Reset cost function.  Used for example if cost function is changed.
function evolve_reset
  integer i
  do i=1,nparams
    cost(i) = evaluate(x1(i,))
  enddo
endf

# --- Do the optimization
function evolve(;gen_max)
  default(gen_max) = 10
  integer i,j,k
  integer a,b,c
  integer count
  chameleon score = cost(1)

  # --- Loop over the generations
  do count=1,gen_max

    # --- loop through population
    do i=1,npop

      # Mutate/Recombine

      # --- Randomly pick three vectors different from each other and 'i'.
      do;a=ranf()*npop+1;until (a ~= i)
      do;b=ranf()*npop+1;until (b ~= i & b ~= a)
      do;c=ranf()*npop+1;until (c ~= i & c ~= a & c ~= b)

      # --- Randomly pick the first parameter
      j = ranf()*nparams+1

      # --- Load parameters into trial, performing binomial trials
      do k=1,nparams
        if (ranf() < crossover | k == nparams) then
          # --- Source for trial is a random vector plus weighted differential
          # --- The last parameter always comes from noisy vector
          trial(j) = x1(c,j) + f*(x1(a,j) - x1(b,j))
        else
          # --- Trial parameter come from x1 itself
          trial(j) = x1(i,j)
        endif
        # --- get next parameter
        j = mod(j,nparams) + 1
      enddo

      # Evaluate/Select

      # --- Evaluate trial function
      score = evaluate(trial)

      if (score <= cost(i)) then
        # --- If trial improves on x1, move trial to secondary array
        # --- and save the new score
        x2(i,) = trial
        cost(i) = score
      else
        # --- otherwise move x1 to secondary array
        x2(i,) = x1(i,)
      endif

    enddo
  
    # --- End of population loop, so copy new parameters into x1
    x1 = x2

  enddo

endf
