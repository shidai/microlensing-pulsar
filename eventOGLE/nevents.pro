; Poisson distribution

p_pulsar =  10^(1.2)*1e-8/1e9*(120e3)
;Calculated assuming OGLE-III NS event rate is 10^(1.2)*1e-8

N_star = 340e6
;OGLE-III monitors 340 million stars

nu = p_pulsar*N_star   

nu = nu*9    ;9 year survey so multiply by 9

;nu = nu*4/9*10   ;OGLE-IV is 4 times greater and 10 years

p = 0
n_tot = 0

FOR n = 1, 20 DO BEGIN

  pn = nu^n*exp(-nu)/factorial(n)
  n_tot = n_tot + n*pn
  p = p + pn

  print, n, pn

ENDFOR

print, 'Probability of a least one event =',  p*100, ' per cent'
print, 'Expected number of events = ', n_tot

END
