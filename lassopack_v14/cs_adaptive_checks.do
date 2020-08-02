prog drop lassoutils
lasso2 $model , l(10) alph(0)
mat b = e(betaAll)
lasso2 $model, l(10) adaptive adal(b)
lasso2 $model, l(10) adaptive adal(b) prestd

lasso2 $model, lcount(2) adaptive adal(b)
lasso2 $model, lcount(2) adaptive adal(b) prestd

lasso2 $model, adaptive adal(b)
lasso2 $model, adaptive adal(b) prestd

lasso2 $model, adaptive
lasso2 $model, adaptive prestd

lasso2 $model, lambda(10) adaptive
lasso2 $model, lambda(10) adaptive prestd

lasso2 $model, lambda(10) adal(b) adaptive
lasso2 $model, lambda(10) adal(b) adaptive prestd

lasso2 $model, lambda(1) adal(b) adatheta(2) adaptive
lasso2 $model, lambda(1) adal(b) adatheta(2) adaptive prestd

lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive prestd
lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive

lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(2) prestd
lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(2)

lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(0.5) prestd
lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(0.5)

lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(3) prestd
lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(2) adaptive adatheta(3)
