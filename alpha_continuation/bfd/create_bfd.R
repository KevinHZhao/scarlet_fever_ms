## this script was made with chatgpt assistance

# read file
lines <- readLines("../../breaks/avg_beta_parms.txt")

# find where each segment starts
seg_idx <- grep("^Segment", lines)

# add end index for easier slicing
seg_idx <- c(seg_idx, length(lines) + 1)

# function to extract parameters from one segment
parse_segment <- function(seg_lines) {
  vals <- seg_lines[grepl("=", seg_lines)]
  vals <- vals[!grepl("^b0", vals)]  # remove b0

  # split into named vector
  parts <- strsplit(vals, "=")
  out <- sapply(parts, function(x) as.numeric(x[2]))
  names(out) <- sapply(parts, function(x) x[1])

  return(out)
}

# template function
make_ode <- function(params, i) {
  sprintf(
    "## For the Forcing Function
# p determines shape of forcing function
# where p = 0 is macpan, p = 1 is sinusoidal
par p=0

s1=%s
s2=%s
s3=%s
c1=%s
c2=%s
c3=%s

par mu=0.02, gamma=24.33, Rzero=7, amp=1

# MACPAN forcing

macpan(t)=s1*sin(2*pi*t)+s2*sin(4*pi*t)+s3*sin(6*pi*t)+c1*cos(2*pi*t)+c2*cos(4*pi*t)+c3*cos(6*pi*t)

# SINUSOIDAL forcing

sinusoid(t)=cos(2*pi*t)

# MACPAN to SINUSOIDAL

transform(t,p)=(1-p)*macpan(t)+p*sinusoid(t)

beta0=Rzero*(gamma+mu)
beta=beta0*(1+amp*transform(t,p))

s'=mu-beta*s*i-mu*s
i'=beta*s*i-(gamma+mu)*i

## INITIAL CONDITIONS:
init S=0.9, I=0.001
## AUXILIARY VARIABLES:
aux R0=Rzero
aux log10s=log10(s)
aux log10i=log10(i)

## PLOT OPTIONS:
@ xp=R0, yp=log10i
@ xlo=0, xhi=40, yhi=0, ylo=-25
@ back=white

## POINCARE MAP SET UP:
@ poimap=section,poivar=t,poipln=1

## range set up
@ range=1, rangeover=Rzero, rangestep=3000
@ rangelow=0, rangehigh=30, rangereset=no

## INTEGRATION OPTIONS:
@ total=650,
@ trans=600
@ dt=0.001

## STORAGE and DATA SAVING OPTIONS
@ maxstor=2000000
@ output=bruteforce_%d.dat
done
",
    params["s1"], params["s2"], params["s3"],
    params["c1"], params["c2"], params["c3"],
    i
  )
}

# loop through segments
for (i in seq_len(length(seg_idx) - 1)) {
  seg_lines <- lines[(seg_idx[i] + 1):(seg_idx[i + 1] - 1)]
  params <- parse_segment(seg_lines)

  ode_text <- make_ode(params, i)

  writeLines(ode_text, sprintf("bruteforce_%d.ode", i))
}
