# pestifer.scripters: 00-01-00_psfgen-build.tcl
####################### Created Fri Apr 17 12:31:09 2026 #######################
package require PestiferCRot
namespace import PestiferCRot::*
package require psfgen
psfcontext mixedcase
topology top_all36_prot.rtf
topology top_all36_cgenff.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_all36_na.rtf
topology toppar_water_ions.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_prot_modify_res.str
topology toppar_all36_moreions.str
topology top_all35_ethers.rtf
pdbalias atom ILE CD1 CD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O
pdbalias atom VCG C01 C1
pdbalias atom VCG C02 C2
pdbalias atom VCG C03 C3
pdbalias atom VCG C04 C4
pdbalias atom VCG C05 C5
pdbalias atom VCG C06 C6
pdbalias atom VCG C07 C7
pdbalias atom VCG C08 C8
pdbalias atom VCG C09 C9
pdbalias atom TIP3 O OH2
pdbalias atom ILE CD1 CD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O
pdbalias atom VCG C01 C1
pdbalias atom VCG C02 C2
pdbalias atom VCG C03 C3
pdbalias atom VCG C04 C4
pdbalias atom VCG C05 C5
pdbalias atom VCG C06 C6
pdbalias atom VCG C07 C7
pdbalias atom VCG C08 C8
pdbalias atom VCG C09 C9
pdbalias atom TIP3 O OH2
pdbalias residue HIS HSD
pdbalias residue PO4 H2PO4
pdbalias residue H2PO H2PO4
pdbalias residue MAN AMAN
pdbalias residue BMA BMAN
pdbalias residue BGLC BGLCNA
pdbalias residue NAG BGLCNA
pdbalias residue NDG AGLCNA
pdbalias residue FUC AFUC
pdbalias residue FUL BFUC
pdbalias residue GAL BGAL
pdbalias residue SIA ANE5AC
pdbalias residue ANE5 ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
pdbalias residue C6DH C6DHPC
pdbalias residue C7DH C7DHPC
pdbalias residue DT THY
pdbalias residue DA ADE
pdbalias residue DC CYT
pdbalias residue DG GUA
pdbalias residue DU URA
pdbalias residue HEM HEME
pdbalias residue TOCL TOCL1
pdbalias residue HIS HSD
pdbalias residue PO4 H2PO4
pdbalias residue H2PO H2PO4
pdbalias residue MAN AMAN
pdbalias residue BMA BMAN
pdbalias residue BGLC BGLCNA
pdbalias residue NAG BGLCNA
pdbalias residue NDG AGLCNA
pdbalias residue FUC AFUC
pdbalias residue FUL BFUC
pdbalias residue GAL BGAL
pdbalias residue SIA ANE5AC
pdbalias residue ANE5 ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
pdbalias residue C6DH C6DHPC
pdbalias residue C7DH C7DHPC
pdbalias residue DT THY
pdbalias residue DA ADE
pdbalias residue DC CYT
pdbalias residue DG GUA
pdbalias residue DU URA
pdbalias residue HEM HEME
pdbalias residue TOCL TOCL1
mol new 4tvp.cif waitfor all
set m1 [molinfo top get id]
set nf [molinfo $m1 get numframes]
if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $m1 }
# Resetting chains and resids for this CIF-source molecule
set TMP [atomselect $m1 "serial 1 to 3565 11345 to 12045"]
$TMP set chain A
$TMP set resid [ list 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 12 12 12 12 12 12 12 13 13 13 13 13 13 13 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 18 18 18 18 18 19 19 19 19 19 19 19 19 19 20 20 20 20 20 20 20 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 23 23 23 23 23 23 23 23 23 23 23 24 24 24 24 24 24 25 25 25 25 25 26 26 26 26 26 26 27 27 27 27 27 27 27 27 28 28 28 28 28 29 29 29 29 29 29 29 29 29 30 30 30 30 30 31 31 31 31 31 31 31 31 31 31 31 31 32 32 32 32 32 32 32 32 32 33 33 33 33 33 33 33 34 34 34 34 34 34 34 34 34 35 35 35 35 35 35 35 35 35 36 36 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 37 38 38 38 38 38 38 38 39 39 39 39 39 39 39 39 39 39 39 39 39 39 40 40 40 40 40 41 41 41 41 41 41 41 42 42 42 42 42 42 42 42 42 42 43 43 43 43 43 44 44 44 44 44 44 45 45 45 45 45 45 45 46 46 46 46 46 46 46 47 47 47 47 47 47 47 48 48 48 48 48 48 48 48 49 49 49 49 49 49 49 50 50 50 50 50 50 50 50 51 51 51 51 51 51 51 52 52 52 52 52 52 52 52 52 53 53 53 53 53 53 53 53 53 54 54 54 54 54 54 54 54 55 55 55 55 55 55 55 55 55 55 56 56 56 56 56 56 56 56 57 57 57 57 57 57 57 57 57 58 58 58 58 58 58 58 58 59 59 59 59 59 59 59 60 60 60 60 60 60 60 61 61 61 61 61 61 61 61 61 62 62 62 62 62 62 62 62 62 63 63 63 63 63 63 63 63 63 63 63 64 64 64 64 64 64 64 64 65 65 65 65 65 65 65 65 66 66 66 66 66 66 66 66 66 66 66 66 66 66 67 67 67 67 67 67 67 67 67 68 68 68 68 68 68 68 68 69 69 69 69 69 69 69 69 70 70 70 70 70 70 70 70 71 71 71 71 71 71 71 72 72 72 72 72 72 72 72 72 73 73 73 73 73 73 73 73 73 74 74 74 74 74 74 74 74 75 75 75 75 75 75 75 75 75 75 76 76 76 76 76 76 76 77 77 77 77 77 77 77 77 78 78 78 78 78 78 78 78 79 79 79 79 79 79 79 79 80 80 80 80 80 80 81 81 81 81 81 81 81 81 82 82 82 82 82 82 82 82 82 82 82 82 82 82 83 83 83 83 83 83 83 83 84 84 84 84 84 84 84 84 84 85 85 85 85 85 85 86 86 86 86 86 86 86 86 87 87 87 87 87 87 87 87 87 88 88 88 88 88 88 88 89 89 89 89 89 89 90 90 90 90 90 90 90 91 91 91 91 91 91 91 91 91 92 92 92 92 92 92 92 92 93 93 93 93 93 93 93 94 94 94 94 94 94 94 95 95 95 95 95 95 95 95 96 96 96 96 96 96 97 97 97 97 97 97 97 98 98 98 98 98 98 98 99 99 99 99 99 99 99 99 100 100 100 100 100 100 100 100 100 101 101 101 101 101 101 102 102 102 102 102 102 102 103 103 103 103 103 103 103 103 104 104 104 104 104 104 104 105 105 105 105 105 105 105 106 106 106 106 106 106 106 106 107 107 107 107 107 107 107 107 108 108 108 108 108 108 108 108 109 109 109 109 109 109 109 110 110 110 110 110 110 110 110 111 111 111 111 111 111 111 111 112 112 112 112 112 112 112 112 113 113 113 113 113 113 113 113 113 113 113 114 114 114 114 115 115 115 115 115 115 115 115 115 116 116 116 116 116 116 116 116 117 117 117 117 117 117 117 117 117 118 118 118 118 118 118 118 118 119 119 119 119 119 119 120 120 120 120 120 120 121 121 121 121 121 121 121 121 121 121 121 122 122 122 122 122 122 122 122 123 123 123 123 123 123 123 123 124 124 124 124 124 124 124 125 125 125 125 125 125 125 126 126 126 126 126 126 126 126 126 127 127 127 127 127 127 127 127 128 128 128 128 128 128 128 128 128 128 128 129 129 129 129 129 129 129 129 130 130 130 130 130 130 130 130 130 131 131 131 131 131 131 131 131 131 132 132 132 132 132 132 132 132 132 133 133 133 133 133 133 133 133 133 134 134 134 134 134 134 134 135 135 135 135 135 135 135 135 135 135 135 135 136 136 136 136 136 136 137 137 137 137 137 137 137 137 138 138 138 138 138 138 138 138 138 138 138 139 139 139 139 139 139 139 139 139 139 139 139 140 140 140 140 140 140 140 140 140 140 140 141 141 141 141 141 141 141 141 142 142 142 142 142 142 142 142 143 143 143 143 143 143 143 144 144 144 144 144 144 144 145 145 145 145 145 145 145 145 145 146 146 146 146 146 146 146 146 147 147 147 147 147 147 147 147 157 157 157 157 157 157 158 158 158 158 158 158 158 158 159 159 159 159 159 159 159 159 159 160 160 160 160 160 160 160 160 160 161 161 161 161 161 161 161 161 161 161 161 161 162 162 162 162 162 162 162 162 162 162 162 163 163 163 163 163 163 163 163 164 164 164 164 164 164 164 164 165 165 165 165 165 165 165 165 166 166 166 166 166 166 167 167 167 167 167 167 167 167 168 168 168 168 168 168 168 169 169 169 169 169 169 170 170 170 170 170 171 171 171 171 171 171 171 171 172 172 172 172 172 172 172 173 173 173 173 173 173 173 173 173 174 174 174 174 174 175 175 175 175 175 175 176 176 176 176 176 176 176 177 177 177 177 177 177 177 177 177 178 178 178 178 178 178 178 179 179 179 179 179 179 180 180 180 180 180 180 180 180 180 180 180 181 181 181 181 181 181 181 181 181 182 182 182 182 182 182 182 183 183 183 183 183 183 183 183 184 184 184 184 184 184 184 185 185 185 185 185 185 185 185 186 186 186 186 186 186 186 186 186 186 187 187 187 187 187 187 187 187 187 187 187 187 188 188 188 188 188 188 189 189 189 189 189 190 190 190 190 190 190 190 191 191 191 191 191 192 192 192 192 193 193 193 193 193 193 193 193 193 193 193 194 194 194 194 194 195 195 195 195 195 195 195 195 196 196 196 196 196 196 196 196 197 197 197 197 197 197 197 197 197 198 198 198 198 198 198 199 199 199 199 199 199 199 199 199 200 200 200 200 200 200 200 200 201 201 201 201 201 201 201 201 201 202 202 202 202 202 202 202 202 202 203 203 203 203 203 203 203 203 203 203 203 204 204 204 204 204 204 204 204 205 205 205 205 206 206 206 206 206 206 206 207 207 207 207 208 208 208 208 208 208 208 209 209 209 209 209 209 210 210 210 210 210 210 210 211 211 211 211 211 211 212 212 212 212 212 212 212 213 213 213 213 213 213 214 214 214 214 214 214 214 215 215 215 215 215 215 215 216 216 216 216 216 216 216 216 216 217 217 217 217 217 217 218 218 218 218 218 218 218 219 219 219 219 219 219 219 219 219 219 220 220 220 220 221 221 221 221 221 221 221 221 222 222 222 222 222 222 222 222 222 223 223 223 223 223 223 223 224 224 224 224 224 224 224 225 225 225 225 225 225 225 226 226 226 226 226 226 227 227 227 227 227 227 227 228 228 228 228 228 228 228 228 228 229 229 229 229 229 229 229 229 230 230 230 230 230 230 230 230 231 231 231 231 231 231 231 231 232 232 232 232 232 232 232 232 233 233 233 233 234 234 234 234 234 234 235 235 235 235 235 235 235 235 236 236 236 236 236 237 237 237 237 237 237 237 237 237 238 238 238 238 238 238 238 238 238 239 239 239 239 239 239 239 239 239 240 240 240 240 240 240 240 241 241 241 241 241 241 241 241 242 242 242 242 242 242 242 242 243 243 243 243 243 243 243 243 243 243 243 244 244 244 244 244 244 245 245 245 245 245 245 245 245 245 246 246 246 246 246 246 246 246 247 247 247 247 247 247 247 247 248 248 248 248 248 248 248 249 249 249 249 249 249 249 249 250 250 250 250 250 250 250 250 251 251 251 251 251 252 252 252 252 252 252 252 252 252 253 253 253 253 253 253 253 253 254 254 254 254 254 254 254 254 255 255 255 255 255 255 255 255 256 256 256 256 256 256 256 257 257 257 257 257 257 257 257 257 258 258 258 258 258 258 258 258 258 258 258 259 259 259 259 259 259 259 259 260 260 260 260 260 260 260 261 261 261 261 261 261 261 262 262 262 262 262 262 262 263 263 263 263 263 263 263 263 263 264 264 264 264 264 264 264 264 265 265 265 265 265 265 265 265 266 266 266 266 266 266 267 267 267 267 267 267 267 268 268 268 268 268 268 268 268 268 268 268 269 269 269 269 269 269 269 270 270 270 270 270 270 270 270 271 271 271 271 271 271 271 271 272 272 272 272 272 272 272 272 273 273 273 273 273 273 273 274 274 274 274 274 274 274 274 274 274 274 275 275 275 275 275 275 275 275 275 276 276 276 276 276 276 277 277 277 277 277 277 277 277 278 278 278 278 278 278 278 278 278 278 278 279 279 279 279 279 279 279 279 280 280 280 280 281 281 281 281 281 281 281 282 282 282 282 283 283 283 283 283 283 283 283 283 284 284 284 284 284 285 285 285 285 285 285 285 285 285 285 285 286 286 286 286 286 286 286 286 286 286 286 286 287 287 287 287 287 288 288 288 288 288 288 288 289 289 289 289 290A 290A 290A 290A 290A 290A 290A 290A 291 291 291 291 291 291 291 291 292 292 292 292 292 292 292 292 293 293 293 293 294 294 294 294 294 294 294 294 295 295 295 295 295 295 295 295 296 296 296 296 296 296 296 296 296 296 296 297 297 297 297 297 297 297 297 297 298 298 298 298 298 299 299 299 299 299 299 299 299 299 299 300 300 300 300 300 300 301 301 301 301 301 301 301 301 302 302 302 302 302 302 302 303 303 303 303 303 303 304 304 304 304 304 304 304 304 304 305 305 305 305 305 306 306 306 306 306 306 306 307 307 307 307 307 307 307 307 307 307 307 307 307 307 308 308 308 308 308 308 308 308 309 309 309 309 309 309 309 309 309 310 310 310 310 310 310 310 311 311 311 311 311 311 311 311 312 312 312 312 313 313 313 313 313 313 313 313 313 314 314 314 314 314 314 314 315 315 315 315 315 315 315 316 316 316 316 316 316 316 316 316 317 317 317 317 317 317 317 317 317 318 318 318 318 318 318 318 318 319 319 319 319 319 319 319 319 319 319 319 320 320 320 320 320 320 320 320 320 321 321 321 321 321 321 321 321 321 321 322 322 322 322 322 322 322 322 322 322 322 323 323 323 323 324 324 324 324 324 324 324 324 325 325 325 325 325 325 325 325 326 326 326 326 326 326 326 327 327 327 327 327 327 327 327 328 328 328 328 328 328 328 328 329 329 329 329 329 329 329 329 329 329 329 330 330 330 330 330 330 330 330 330 330 330 331 331 331 331 331 332 332 332 332 332 332 332 332 333 333 333 333 333 333 334 334 334 334 334 334 335 335 335 335 336 336 336 336 337 337 337 337 337 337 337 337 338 338 338 338 338 338 338 338 339 339 339 339 339 339 339 339 339 340 340 340 340 340 340 340 341 341 341 341 341 341 341 342 342 342 342 342 342 342 343 343 343 343 343 343 343 343 343 343 344 344 344 344 344 344 345 345 345 345 345 345 345 345 345 345 345 346 346 346 346 346 346 346 346 347 347 347 347 347 347 348 348 348 348 349 349 349 349 350 350 350 350 350 350 350 350 350 351 351 351 351 351 351 351 351 351 351 351 352 352 352 352 352 352 352 352 352 352 352 353 353 353 353 353 353 353 353 353 353 353 353 354 354 354 354 354 354 355 355 355 355 355 355 355 355 356 356 356 356 356 356 356 357 357 357 357 357 357 358 358 358 358 359 359 359 359 359 359 359 359 360 360 360 360 360 360 360 360 360 360 360 361 361 361 361 361 361 361 361 362 362 362 362 362 362 363 363 363 363 363 363 363 364 364 364 364 364 364 364 364 364 364 364 364 364 364 365 365 365 365 365 365 365 365 366 366 366 366 366 366 367 367 367 367 367 367 367 367 379 379 379 379 379 379 379 379 380 380 380 380 380 380 380 380 381 381 381 381 381 381 382 382 382 382 382 382 382 382 383 383 383 383 383 383 383 384 384 384 384 384 384 384 384 385 385 385 385 385 385 385 386 386 386 386 386 386 387 387 387 387 387 387 387 387 387 387 387 388 388 388 388 388 388 388 388 389 389 389 389 389 389 389 389 389 390 390 390 390 390 390 390 390 390 391 391 391 391 391 391 391 391 392 392 392 392 392 392 392 392 393 393 393 393 393 393 393 393 394 394 394 394 394 394 394 394 395 395 395 395 395 395 395 395 395 395 395 395 395 395 396 396 396 396 396 396 396 396 396 397 397 397 397 397 397 397 397 397 397 397 398 398 398 398 398 398 398 398 399 399 399 399 400 400 400 400 400 400 400 400 400 401 401 401 401 401 402 402 402 402 402 402 402 402 403 403 403 403 403 403 403 403 403 403 403 403 404 404 404 404 404 405 405 405 405 405 405 405 406 406 406 406 406 406 406 407 407 407 407 407 407 407 407 408 408 408 408 408 408 408 408 408 409 409 409 409 410 410 410 410 410 410 410 411 411 411 411 411 411 411 411 412 412 412 412 412 412 412 412 412 412 412 413 413 413 413 413 413 414 414 414 414 414 414 414 415 415 415 415 415 415 416 416 416 416 416 416 416 416 417 417 417 417 417 417 417 417 418 418 418 418 418 418 418 419 419 419 419 420 420 420 420 420 420 420 420 421 421 421 421 421 421 421 421 422 422 422 422 422 422 422 422 423 423 423 423 423 423 423 424 424 424 424 424 424 424 424 424 424 424 425 425 425 425 425 425 425 425 426 426 426 426 427 427 427 427 428 428 428 428 428 428 429 429 429 429 429 429 429 430 430 430 430 430 430 430 430 431 431 431 431 431 431 432 432 432 432 432 432 432 433 433 433 433 433 433 433 434 434 434 434 434 434 434 434 434 435 435 435 435 435 435 435 436 436 436 436 436 436 436 436 436 436 436 437 437 437 437 437 437 437 437 437 437 437 438 438 438 438 438 438 438 439 439 439 439 440 440 440 440 441 441 441 441 442 442 442 442 442 442 442 442 443 443 443 443 443 443 443 443 444 444 444 444 444 444 444 444 444 444 444 445 445 445 445 445 445 445 445 446 446 446 446 446 446 446 446 447 447 447 447 447 447 447 447 447 447 447 447 447 447 448 448 448 448 448 448 448 448 448 448 448 449 449 449 449 449 449 450 450 450 450 450 450 450 450 450 451 451 451 451 451 451 451 451 452 452 452 452 452 452 452 452 452 452 452 452 453 453 453 453 453 453 453 453 453 454 454 454 454 454 454 454 454 454 454 454 454 455 455 455 455 455 455 455 455 455 456 456 456 456 456 456 456 457 457 457 457 457 457 457 458 458 458 458 458 458 458 458 458 459 459 459 459 459 459 459 459 460 460 460 460 460 460 460 460 460 461 461 461 461 461 461 461 462 462 462 462 462 462 462 462 463 463 463 463 464 464 464 464 464 464 464 465 465 465 465 465 466 466 466 466 466 466 466 467 467 467 467 467 467 467 468 468 468 468 468 468 468 468 468 468 468 469 469 469 469 469 469 470 470 470 470 470 470 470 470 470 471 471 471 471 471 471 471 471 471 471 471 472 472 472 472 472 472 472 472 472 472 472 473 473 473 473 473 473 473 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 ]
$TMP delete
set TMP [atomselect $m1 "serial 3566 to 4566 12081 to 12122"]
$TMP set chain B
$TMP set resid [ list 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 12 13 13 13 13 14 14 14 14 14 15 15 15 15 15 16 16 16 16 17 17 17 17 17 17 18 18 18 18 18 18 18 19 19 19 19 19 19 19 19 20 20 20 20 21 21 21 21 21 22 22 22 22 22 23 23 23 23 23 23 24 24 24 24 24 24 24 24 25 25 25 25 25 25 25 26 26 26 26 26 26 26 26 27 27 27 27 27 27 27 28 28 28 28 28 28 28 29 29 29 29 29 29 29 29 29 30 30 30 30 30 31 31 31 31 31 31 31 31 31 31 31 32 32 32 32 32 32 32 32 33 33 33 33 33 33 33 33 34 34 34 34 34 34 34 34 35 35 35 35 35 35 36 36 36 36 58 58 58 58 58 58 58 59 59 59 59 59 59 59 60 60 60 60 60 60 60 60 60 60 60 60 60 60 61 61 61 61 62 62 62 62 62 62 62 62 63 63 63 63 63 63 63 63 63 64 64 64 64 64 64 64 64 64 65 65 65 65 65 65 65 65 66 66 66 66 66 66 66 66 66 67 67 67 67 67 68 68 68 68 68 68 68 68 68 68 68 69 69 69 69 69 69 69 70 70 70 70 70 70 70 70 71 71 71 71 71 72 72 72 72 72 72 72 73 73 73 73 73 73 73 73 73 74 74 74 74 74 74 74 74 74 74 74 75 75 75 75 75 75 75 75 75 75 75 75 76 76 76 76 76 76 76 76 77 77 77 77 77 77 77 77 77 77 77 78 78 78 78 78 78 78 78 79 79 79 79 79 79 79 79 79 80 80 80 80 80 80 80 80 80 81 81 81 81 81 81 81 81 82 82 82 82 82 82 82 82 83 83 83 83 84 84 84 84 84 84 84 84 85 85 85 85 85 85 85 85 85 85 85 85 85 85 86 86 86 86 87 87 87 87 87 87 88 88 88 88 88 88 89 89 89 89 90 90 90 90 90 90 90 90 90 91 91 91 91 91 91 91 91 92 92 92 92 92 92 92 92 93 93 93 93 93 93 94 94 94 94 94 94 95 95 95 95 95 95 95 96 96 96 96 96 96 96 96 97 97 97 97 97 97 97 98 98 98 98 98 98 98 99 99 99 99 99 99 99 99 99 99 99 99 99 99 100 100 100 100 100 100 100 100 101 101 101 101 101 101 102 102 102 102 102 102 103 103 103 103 103 103 103 103 103 103 103 103 103 103 104 104 104 104 104 104 105 105 105 105 105 105 105 105 106 106 106 106 106 106 106 106 106 106 106 107 107 107 107 107 107 107 107 108 108 108 108 108 108 108 108 109 109 109 109 109 109 110 110 110 110 110 110 110 110 110 111 111 111 111 111 111 111 111 112 112 112 112 112 112 112 112 112 112 112 112 112 112 113 113 113 113 113 113 113 113 114 114 114 114 114 114 114 114 115 115 115 115 115 115 115 115 116 116 116 116 116 116 116 117 117 117 117 117 117 117 117 117 117 117 117 117 117 118 118 118 118 118 118 118 118 119 119 119 119 119 119 119 119 119 120 120 120 120 120 120 120 120 120 120 120 120 120 120 121 121 121 121 121 121 121 121 122 122 122 122 122 122 122 122 122 123 123 123 123 123 123 123 123 123 124 124 124 124 124 124 124 124 125 125 125 125 125 125 126 126 126 126 126 126 126 126 127 127 127 127 127 127 127 127 127 127 127 127 128 128 128 128 128 128 128 129 129 129 129 129 129 129 129 129 130 130 130 130 130 130 130 130 131 131 131 131 131 131 131 131 132 132 132 132 132 132 132 132 132 132 132 132 133 133 133 133 134 134 134 134 134 134 134 134 135 135 135 135 135 135 135 135 136 136 136 136 136 136 136 136 136 137 137 137 137 137 137 137 137 137 138 138 138 138 138 138 139 139 139 139 139 139 139 139 139 140 140 140 140 140 140 140 140 141 141 141 141 141 141 141 141 141 142 142 142 142 142 142 142 142 142 143 143 143 143 143 143 143 143 143 144 144 144 144 144 144 144 144 144 145 145 145 145 145 145 145 145 146 146 146 146 146 146 146 146 146 147 147 147 147 147 147 147 147 147 148 148 148 148 148 148 148 148 149 149 149 149 149 149 149 149 150 150 150 150 150 150 150 150 151 151 151 151 151 152 152 152 152 152 152 152 152 153 153 153 153 153 153 153 153 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 ]
$TMP delete
set TMP [atomselect $m1 "serial 12046 to 12050"]
$TMP set chain Y
$TMP set resid [ list 601 601 601 601 601 ]
$TMP delete
set TMP [atomselect $m1 "serial 12051 to 12055"]
$TMP set chain Z
$TMP set resid [ list 602 602 602 602 602 ]
$TMP delete
set TMP [atomselect $m1 "serial 12056 to 12060"]
$TMP set chain AA
$TMP set resid [ list 603 603 603 603 603 ]
$TMP delete
set TMP [atomselect $m1 "serial 12061 to 12065"]
$TMP set chain BA
$TMP set resid [ list 604 604 604 604 604 ]
$TMP delete
set TMP [atomselect $m1 "serial 12066 to 12070"]
$TMP set chain CA
$TMP set resid [ list 605 605 605 605 605 ]
$TMP delete
set TMP [atomselect $m1 "serial 12071 to 12075"]
$TMP set chain DA
$TMP set resid [ list 606 606 606 606 606 ]
$TMP delete
set TMP [atomselect $m1 "serial 12076 to 12080"]
$TMP set chain EA
$TMP set resid [ list 607 607 607 607 607 ]
$TMP delete
set TMP [atomselect $m1 "serial 12123 to 12127"]
$TMP set chain IA
$TMP set resid [ list 701 701 701 701 701 ]
$TMP delete
set TMP [atomselect $m1 "serial 12128 to 12132"]
$TMP set chain JA
$TMP set resid [ list 702 702 702 702 702 ]
$TMP delete
set TMP [atomselect $m1 "serial 12152 to 12168"]
$TMP set chain MA
$TMP set resid [ list 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717 ]
$TMP delete
set TMP [atomselect $m1 "serial 12169 12170"]
$TMP set chain NA
$TMP set resid [ list 801 802 ]
$TMP delete
# Done.
mol top $m1
############################## Transform 0 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
# A.U. chain A: Image chain A
# A.U. chain B: Image chain B
# A.U. chain AG01: Image chain AG01
# A.U. chain AG02: Image chain AG02
# A.U. chain AG03: Image chain AG03
# A.U. chain AG04: Image chain AG04
# A.U. chain AG05: Image chain AG05
# A.U. chain AG06: Image chain AG06
# A.U. chain AG07: Image chain AG07
# A.U. chain AG08: Image chain AG08
# A.U. chain AG09: Image chain AG09
# A.U. chain AG10: Image chain AG10
# A.U. chain AG11: Image chain AG11
# A.U. chain AG12: Image chain AG12
# A.U. chain AG13: Image chain AG13
# A.U. chain AG14: Image chain AG14
# A.U. chain AG15: Image chain AG15
# A.U. chain AG16: Image chain AG16
# A.U. chain AG17: Image chain AG17
# A.U. chain AG18: Image chain AG18
# A.U. chain BG01: Image chain BG01
# A.U. chain BG02: Image chain BG02
# A.U. chain BG03: Image chain BG03
# A.U. chain Y: Image chain Y
# A.U. chain Z: Image chain Z
# A.U. chain AA: Image chain AA
# A.U. chain BA: Image chain BA
# A.U. chain CA: Image chain CA
# A.U. chain DA: Image chain DA
# A.U. chain EA: Image chain EA
# A.U. chain IA: Image chain IA
# A.U. chain JA: Image chain JA
# A.U. chain MA: Image chain MA
# A.U. chain NA: Image chain NA
############################### Segments follow ################################
# Segment A begins
############################### Segment A begins ###############################
set A00 [atomselect $m1 "serial 1 to 1196"]
$A00 set segname A
$A00 writepdb segtype_polymer_A_1_to_147.pdb
set A02 [atomselect $m1 "serial 1197 to 2817"]
$A02 set segname A
$A02 writepdb segtype_polymer_A_157_to_367.pdb
set A04 [atomselect $m1 "serial 2818 to 3565"]
$A04 set segname A
$A04 writepdb segtype_polymer_A_379_to_473.pdb
segment A {
    pdb segtype_polymer_A_1_to_147.pdb
    residue 148 GLU A
    residue 149 ASN A
    residue 150 GLN A
    residue 151 GLY A
    residue 152 ASN A
    residue 153 ARG A
    residue 154 SER A
    residue 155 ASN A
    residue 156 ASN A
    residue 156A GLY A
    pdb segtype_polymer_A_157_to_367.pdb
    residue 368 THR A
    residue 369 SER A
    residue 370 VAL A
    residue 371 GLN A
    residue 372 GLY A
    residue 373 SER A
    residue 374 ASN A
    residue 375 SER A
    residue 376 THR A
    residue 377 GLY A
    residue 378 SER A
    residue 378A GLY A
    pdb segtype_polymer_A_379_to_473.pdb
    mutate 469 ALA
}
################################ End segment A #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_A_1_to_147.pdb A
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_A_157_to_367.pdb A
# Subsegment [4] is a resolved run
coordpdb segtype_polymer_A_379_to_473.pdb A
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at A_GLU148 from A_ASN147
coord A 148 N {-269.51259 -129.47827 -46.16939}
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at A_THR368 from A_ASN367
coord A 368 N {-309.98814 -112.22876 -7.07846}
####################### Intra-segmental terminal patches #######################
patch CTER A:156
patch NTER A:157
delatom A 156A
patch CTER A:378
patch NTER A:379
delatom A 378A
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment A ends ################################
# Segment A ends
# Segment B begins
############################### Segment B begins ###############################
set B01 [atomselect $m1 "serial 3566 to 3769"]
$B01 set segname B
$B01 writepdb segtype_polymer_B_7_to_36.pdb
set B03 [atomselect $m1 "serial 3770 to 4566"]
$B03 set segname B
$B03 writepdb segtype_polymer_B_58_to_153.pdb
segment B {
    pdb segtype_polymer_B_7_to_36.pdb
    residue 37 ILE B
    residue 38 VAL B
    residue 39 GLN B
    residue 40 GLN B
    residue 41 GLN B
    residue 42 SER B
    residue 43 ASN B
    residue 44 LEU B
    residue 45 LEU B
    residue 46 ARG B
    residue 47 ALA B
    residue 48 PRO B
    residue 49 GLU B
    residue 50 ALA B
    residue 51 GLN B
    residue 52 GLN B
    residue 53 HSD B
    residue 54 LEU B
    residue 55 LEU B
    residue 56 LYS B
    residue 57 LEU B
    residue 57A GLY B
    pdb segtype_polymer_B_58_to_153.pdb
    mutate 94 THR
    mutate 48 ILE
}
################################ End segment B #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_B_7_to_36.pdb B
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_B_58_to_153.pdb B
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B_ILE37 from B_GLY36
coord B 37 N {-266.96525 -151.71222 18.97426}
####################### Intra-segmental terminal patches #######################
patch CTER B:57
patch NTER B:58
delatom B 57A
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment B ends ################################
# Segment B ends
# Segment AG01 begins
##################### Segment AG01 begins as image of AG01 #####################
set AG01 [atomselect $m1 "serial 11345 to 11427"]
$AG01 set segname AG01
$AG01 set chain A
$AG01 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7]
$AG01 writepdb segtype_generic_AG01.pdb
segment AG01 {
    first none
    last none
    pdb segtype_generic_AG01.pdb
}
coordpdb segtype_generic_AG01.pdb AG01
############################## Segment AG01 ends ###############################
# Segment AG01 ends
# Segment AG02 begins
##################### Segment AG02 begins as image of AG02 #####################
set AG02 [atomselect $m1 "serial 11428 to 11477"]
$AG02 set segname AG02
$AG02 set chain A
$AG02 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4]
$AG02 writepdb segtype_generic_AG02.pdb
segment AG02 {
    first none
    last none
    pdb segtype_generic_AG02.pdb
}
coordpdb segtype_generic_AG02.pdb AG02
############################## Segment AG02 ends ###############################
# Segment AG02 ends
# Segment AG03 begins
##################### Segment AG03 begins as image of AG03 #####################
set AG03 [atomselect $m1 "serial 11478 to 11505"]
$AG03 set segname AG03
$AG03 set chain A
$AG03 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG03 writepdb segtype_generic_AG03.pdb
segment AG03 {
    first none
    last none
    pdb segtype_generic_AG03.pdb
}
coordpdb segtype_generic_AG03.pdb AG03
############################## Segment AG03 ends ###############################
# Segment AG03 ends
# Segment AG04 begins
##################### Segment AG04 begins as image of AG04 #####################
set AG04 [atomselect $m1 "serial 11506 to 11577"]
$AG04 set segname AG04
$AG04 set chain A
$AG04 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6]
$AG04 writepdb segtype_generic_AG04.pdb
segment AG04 {
    first none
    last none
    pdb segtype_generic_AG04.pdb
}
coordpdb segtype_generic_AG04.pdb AG04
############################## Segment AG04 ends ###############################
# Segment AG04 ends
# Segment AG05 begins
##################### Segment AG05 begins as image of AG05 #####################
set AG05 [atomselect $m1 "serial 11578 to 11605"]
$AG05 set segname AG05
$AG05 set chain A
$AG05 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG05 writepdb segtype_generic_AG05.pdb
segment AG05 {
    first none
    last none
    pdb segtype_generic_AG05.pdb
}
coordpdb segtype_generic_AG05.pdb AG05
############################## Segment AG05 ends ###############################
# Segment AG05 ends
# Segment AG06 begins
##################### Segment AG06 begins as image of AG06 #####################
set AG06 [atomselect $m1 "serial 11606 to 11633"]
$AG06 set segname AG06
$AG06 set chain A
$AG06 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG06 writepdb segtype_generic_AG06.pdb
segment AG06 {
    first none
    last none
    pdb segtype_generic_AG06.pdb
}
coordpdb segtype_generic_AG06.pdb AG06
############################## Segment AG06 ends ###############################
# Segment AG06 ends
# Segment AG07 begins
##################### Segment AG07 begins as image of AG07 #####################
set AG07 [atomselect $m1 "serial 11634 to 11661"]
$AG07 set segname AG07
$AG07 set chain A
$AG07 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG07 writepdb segtype_generic_AG07.pdb
segment AG07 {
    first none
    last none
    pdb segtype_generic_AG07.pdb
}
coordpdb segtype_generic_AG07.pdb AG07
############################## Segment AG07 ends ###############################
# Segment AG07 ends
# Segment AG08 begins
##################### Segment AG08 begins as image of AG08 #####################
set AG08 [atomselect $m1 "serial 11662 to 11700"]
$AG08 set segname AG08
$AG08 set chain A
$AG08 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3]
$AG08 writepdb segtype_generic_AG08.pdb
segment AG08 {
    first none
    last none
    pdb segtype_generic_AG08.pdb
}
coordpdb segtype_generic_AG08.pdb AG08
############################## Segment AG08 ends ###############################
# Segment AG08 ends
# Segment AG09 begins
##################### Segment AG09 begins as image of AG09 #####################
set AG09 [atomselect $m1 "serial 11701 to 11761"]
$AG09 set segname AG09
$AG09 set chain A
$AG09 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5]
$AG09 writepdb segtype_generic_AG09.pdb
segment AG09 {
    first none
    last none
    pdb segtype_generic_AG09.pdb
}
coordpdb segtype_generic_AG09.pdb AG09
############################## Segment AG09 ends ###############################
# Segment AG09 ends
# Segment AG10 begins
##################### Segment AG10 begins as image of AG10 #####################
set AG10 [atomselect $m1 "serial 11762 to 11789"]
$AG10 set segname AG10
$AG10 set chain A
$AG10 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG10 writepdb segtype_generic_AG10.pdb
segment AG10 {
    first none
    last none
    pdb segtype_generic_AG10.pdb
}
coordpdb segtype_generic_AG10.pdb AG10
############################## Segment AG10 ends ###############################
# Segment AG10 ends
# Segment AG11 begins
##################### Segment AG11 begins as image of AG11 #####################
set AG11 [atomselect $m1 "serial 11790 to 11817"]
$AG11 set segname AG11
$AG11 set chain A
$AG11 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG11 writepdb segtype_generic_AG11.pdb
segment AG11 {
    first none
    last none
    pdb segtype_generic_AG11.pdb
}
coordpdb segtype_generic_AG11.pdb AG11
############################## Segment AG11 ends ###############################
# Segment AG11 ends
# Segment AG12 begins
##################### Segment AG12 begins as image of AG12 #####################
set AG12 [atomselect $m1 "serial 11818 to 11845"]
$AG12 set segname AG12
$AG12 set chain A
$AG12 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG12 writepdb segtype_generic_AG12.pdb
segment AG12 {
    first none
    last none
    pdb segtype_generic_AG12.pdb
}
coordpdb segtype_generic_AG12.pdb AG12
############################## Segment AG12 ends ###############################
# Segment AG12 ends
# Segment AG13 begins
##################### Segment AG13 begins as image of AG13 #####################
set AG13 [atomselect $m1 "serial 11846 to 11873"]
$AG13 set segname AG13
$AG13 set chain A
$AG13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$AG13 writepdb segtype_generic_AG13.pdb
segment AG13 {
    first none
    last none
    pdb segtype_generic_AG13.pdb
}
coordpdb segtype_generic_AG13.pdb AG13
############################## Segment AG13 ends ###############################
# Segment AG13 ends
# Segment AG14 begins
##################### Segment AG14 begins as image of AG14 #####################
set AG14 [atomselect $m1 "serial 11874 to 11989"]
$AG14 set segname AG14
$AG14 set chain A
$AG14 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10]
$AG14 writepdb segtype_generic_AG14.pdb
segment AG14 {
    first none
    last none
    pdb segtype_generic_AG14.pdb
}
coordpdb segtype_generic_AG14.pdb AG14
############################## Segment AG14 ends ###############################
# Segment AG14 ends
# Segment AG15 begins
##################### Segment AG15 begins as image of AG15 #####################
set AG15 [atomselect $m1 "serial 11990 to 12003"]
$AG15 set segname AG15
$AG15 set chain A
$AG15 set resid [list 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276]
$AG15 writepdb segtype_generic_AG15.pdb
segment AG15 {
    first none
    last none
    pdb segtype_generic_AG15.pdb
}
coordpdb segtype_generic_AG15.pdb AG15
############################## Segment AG15 ends ###############################
# Segment AG15 ends
# Segment AG16 begins
##################### Segment AG16 begins as image of AG16 #####################
set AG16 [atomselect $m1 "serial 12004 to 12017"]
$AG16 set segname AG16
$AG16 set chain A
$AG16 set resid [list 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839]
$AG16 writepdb segtype_generic_AG16.pdb
segment AG16 {
    first none
    last none
    pdb segtype_generic_AG16.pdb
}
coordpdb segtype_generic_AG16.pdb AG16
############################## Segment AG16 ends ###############################
# Segment AG16 ends
# Segment AG17 begins
##################### Segment AG17 begins as image of AG17 #####################
set AG17 [atomselect $m1 "serial 12018 to 12031"]
$AG17 set segname AG17
$AG17 set chain A
$AG17 set resid [list 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355]
$AG17 writepdb segtype_generic_AG17.pdb
segment AG17 {
    first none
    last none
    pdb segtype_generic_AG17.pdb
}
coordpdb segtype_generic_AG17.pdb AG17
############################## Segment AG17 ends ###############################
# Segment AG17 ends
# Segment AG18 begins
##################### Segment AG18 begins as image of AG18 #####################
set AG18 [atomselect $m1 "serial 12032 to 12045"]
$AG18 set segname AG18
$AG18 set chain A
$AG18 set resid [list 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133]
$AG18 writepdb segtype_generic_AG18.pdb
segment AG18 {
    first none
    last none
    pdb segtype_generic_AG18.pdb
}
coordpdb segtype_generic_AG18.pdb AG18
############################## Segment AG18 ends ###############################
# Segment AG18 ends
# Segment BG01 begins
##################### Segment BG01 begins as image of BG01 #####################
set BG01 [atomselect $m1 "serial 12081 to 12094"]
$BG01 set segname BG01
$BG01 set chain B
$BG01 set resid [list 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611]
$BG01 writepdb segtype_generic_BG01.pdb
segment BG01 {
    first none
    last none
    pdb segtype_generic_BG01.pdb
}
coordpdb segtype_generic_BG01.pdb BG01
############################## Segment BG01 ends ###############################
# Segment BG01 ends
# Segment BG02 begins
##################### Segment BG02 begins as image of BG02 #####################
set BG02 [atomselect $m1 "serial 12095 to 12108"]
$BG02 set segname BG02
$BG02 set chain B
$BG02 set resid [list 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618]
$BG02 writepdb segtype_generic_BG02.pdb
segment BG02 {
    first none
    last none
    pdb segtype_generic_BG02.pdb
}
coordpdb segtype_generic_BG02.pdb BG02
############################## Segment BG02 ends ###############################
# Segment BG02 ends
# Segment BG03 begins
##################### Segment BG03 begins as image of BG03 #####################
set BG03 [atomselect $m1 "serial 12109 to 12122"]
$BG03 set segname BG03
$BG03 set chain B
$BG03 set resid [list 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637]
$BG03 writepdb segtype_generic_BG03.pdb
segment BG03 {
    first none
    last none
    pdb segtype_generic_BG03.pdb
}
coordpdb segtype_generic_BG03.pdb BG03
############################## Segment BG03 ends ###############################
# Segment BG03 ends
# Segment Y begins
######################## Segment Y begins as image of Y ########################
set Y [atomselect $m1 "serial 12046 to 12050"]
$Y set segname Y
$Y set resid [list 601 601 601 601 601]
$Y writepdb segtype_generic_Y.pdb
segment Y {
    first none
    last none
    pdb segtype_generic_Y.pdb
}
coordpdb segtype_generic_Y.pdb Y
################################ Segment Y ends ################################
# Segment Y ends
# Segment Z begins
######################## Segment Z begins as image of Z ########################
set Z [atomselect $m1 "serial 12051 to 12055"]
$Z set segname Z
$Z set resid [list 602 602 602 602 602]
$Z writepdb segtype_generic_Z.pdb
segment Z {
    first none
    last none
    pdb segtype_generic_Z.pdb
}
coordpdb segtype_generic_Z.pdb Z
################################ Segment Z ends ################################
# Segment Z ends
# Segment AA begins
####################### Segment AA begins as image of AA #######################
set AA [atomselect $m1 "serial 12056 to 12060"]
$AA set segname AA
$AA set resid [list 603 603 603 603 603]
$AA writepdb segtype_generic_AA.pdb
segment AA {
    first none
    last none
    pdb segtype_generic_AA.pdb
}
coordpdb segtype_generic_AA.pdb AA
############################### Segment AA ends ################################
# Segment AA ends
# Segment BA begins
####################### Segment BA begins as image of BA #######################
set BA [atomselect $m1 "serial 12061 to 12065"]
$BA set segname BA
$BA set resid [list 604 604 604 604 604]
$BA writepdb segtype_generic_BA.pdb
segment BA {
    first none
    last none
    pdb segtype_generic_BA.pdb
}
coordpdb segtype_generic_BA.pdb BA
############################### Segment BA ends ################################
# Segment BA ends
# Segment CA begins
####################### Segment CA begins as image of CA #######################
set CA [atomselect $m1 "serial 12066 to 12070"]
$CA set segname CA
$CA set resid [list 605 605 605 605 605]
$CA writepdb segtype_generic_CA.pdb
segment CA {
    first none
    last none
    pdb segtype_generic_CA.pdb
}
coordpdb segtype_generic_CA.pdb CA
############################### Segment CA ends ################################
# Segment CA ends
# Segment DA begins
####################### Segment DA begins as image of DA #######################
set DA [atomselect $m1 "serial 12071 to 12075"]
$DA set segname DA
$DA set resid [list 606 606 606 606 606]
$DA writepdb segtype_generic_DA.pdb
segment DA {
    first none
    last none
    pdb segtype_generic_DA.pdb
}
coordpdb segtype_generic_DA.pdb DA
############################### Segment DA ends ################################
# Segment DA ends
# Segment EA begins
####################### Segment EA begins as image of EA #######################
set EA [atomselect $m1 "serial 12076 to 12080"]
$EA set segname EA
$EA set resid [list 607 607 607 607 607]
$EA writepdb segtype_generic_EA.pdb
segment EA {
    first none
    last none
    pdb segtype_generic_EA.pdb
}
coordpdb segtype_generic_EA.pdb EA
############################### Segment EA ends ################################
# Segment EA ends
# Segment IA begins
####################### Segment IA begins as image of IA #######################
set IA [atomselect $m1 "serial 12123 to 12127"]
$IA set segname IA
$IA set resid [list 701 701 701 701 701]
$IA writepdb segtype_generic_IA.pdb
segment IA {
    first none
    last none
    pdb segtype_generic_IA.pdb
}
coordpdb segtype_generic_IA.pdb IA
############################### Segment IA ends ################################
# Segment IA ends
# Segment JA begins
####################### Segment JA begins as image of JA #######################
set JA [atomselect $m1 "serial 12128 to 12132"]
$JA set segname JA
$JA set resid [list 702 702 702 702 702]
$JA writepdb segtype_generic_JA.pdb
segment JA {
    first none
    last none
    pdb segtype_generic_JA.pdb
}
coordpdb segtype_generic_JA.pdb JA
############################### Segment JA ends ################################
# Segment JA ends
# Segment MA begins
####################### Segment MA begins as image of MA #######################
set MA [atomselect $m1 "serial 12152 to 12168"]
$MA set segname MA
$MA set resid [list 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717]
$MA writepdb segtype_generic_MA.pdb
segment MA {
    first none
    last none
    pdb segtype_generic_MA.pdb
}
coordpdb segtype_generic_MA.pdb MA
############################### Segment MA ends ################################
# Segment MA ends
# Segment NA begins
####################### Segment NA begins as image of NA #######################
set NA [atomselect $m1 "serial 12169 12170"]
$NA set segname NA
$NA set resid [list 801 802]
$NA writepdb segtype_generic_NA.pdb
segment NA {
    first none
    last none
    pdb segtype_generic_NA.pdb
}
coordpdb segtype_generic_NA.pdb NA
############################### Segment NA ends ################################
# Segment NA ends
############################# DISU patches follow ##############################
patch DISU A:24 A:44
patch DISU A:89 A:175
patch DISU A:96 A:166
patch DISU A:101 A:119
patch DISU A:188 A:217
patch DISU A:198 A:209
patch DISU A:266 A:300
patch DISU A:347 A:413
patch DISU A:354 A:386
patch DISU B:87 B:93
############################# LINK patches follow ##############################
patch NGLA A:58 AG01:1
patch NGLB A:103 AG18:1133
patch NGLA A:107 AG02:1
patch NGLA A:118 AG09:1
patch NGLA A:122 AG10:1
patch NGLA A:167 AG11:1
patch NGLA A:204 AG03:1
patch NGLA A:232 AG04:1
patch NGLA A:246 AG15:1276
patch NGLA A:265 AG12:1
patch NGLB A:271 AG13:1
patch NGLA A:301 AG14:1
patch NGLB A:308 AG16:1839
patch NGLA A:324 AG17:1355
patch NGLA A:332 AG05:1
patch NGLA A:355 AG06:1
patch NGLB A:361 AG07:1
patch NGLA A:416 AG08:1
patch NGLB B:100 BG01:1611
patch NGLA B:107 BG02:1618
patch NGLA B:126 BG03:1637
patch 14ab AG01:1 AG01:2
patch 14ab AG01:2 AG01:3
patch 16BT AG01:3 AG01:4
patch 13aa AG01:3 AG01:7
patch 13ba AG01:4 AG01:5
patch 16AT AG01:4 AG01:6
patch 14ab AG02:1 AG02:2
patch 14ab AG02:2 AG02:3
patch 13ba AG02:3 AG02:4
patch 14ab AG03:1 AG03:2
patch 14bb AG04:1 AG04:2
patch 14aa AG04:2 AG04:3
patch 13ba AG04:3 AG04:4
patch 16AT AG04:3 AG04:6
patch 12ba AG04:4 AG04:5
patch 14bb AG05:1 AG05:2
patch 14bb AG06:1 AG06:2
patch 14bb AG07:1 AG07:2
patch 14ab AG08:1 AG08:2
patch 14bb AG08:2 AG08:3
patch 14ab AG09:1 AG09:2
patch 14ab AG09:2 AG09:3
patch 13ba AG09:3 AG09:4
patch 16AT AG09:3 AG09:5
patch 14ab AG10:1 AG10:2
patch 14ab AG11:1 AG11:2
patch 14ab AG12:1 AG12:2
patch 14bb AG13:1 AG13:2
patch 14bb AG14:1 AG14:2
patch 14bb AG14:2 AG14:3
patch 13aa AG14:3 AG14:4
patch 16AT AG14:3 AG14:7
patch 12aa AG14:4 AG14:5
patch 12ab AG14:5 AG14:6
patch 13bb AG14:7 AG14:8
patch 16AT AG14:7 AG14:10
patch 12bb AG14:8 AG14:9
############################### Transform 0 ends ###############################
############################## Transform 1 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
# A.U. chain A: Image chain C
# A.U. chain B: Image chain D
# A.U. chain AG01: Image chain CG01
# A.U. chain AG02: Image chain CG02
# A.U. chain AG03: Image chain CG03
# A.U. chain AG04: Image chain CG04
# A.U. chain AG05: Image chain CG05
# A.U. chain AG06: Image chain CG06
# A.U. chain AG07: Image chain CG07
# A.U. chain AG08: Image chain CG08
# A.U. chain AG09: Image chain CG09
# A.U. chain AG10: Image chain CG10
# A.U. chain AG11: Image chain CG11
# A.U. chain AG12: Image chain CG12
# A.U. chain AG13: Image chain CG13
# A.U. chain AG14: Image chain CG14
# A.U. chain AG15: Image chain CG15
# A.U. chain AG16: Image chain CG16
# A.U. chain AG17: Image chain CG17
# A.U. chain AG18: Image chain CG18
# A.U. chain BG01: Image chain DG01
# A.U. chain BG02: Image chain DG02
# A.U. chain BG03: Image chain DG03
# A.U. chain Y: Image chain HB
# A.U. chain Z: Image chain IB
# A.U. chain AA: Image chain JB
# A.U. chain BA: Image chain KB
# A.U. chain CA: Image chain LB
# A.U. chain DA: Image chain MB
# A.U. chain EA: Image chain NB
# A.U. chain IA: Image chain OB
# A.U. chain JA: Image chain PB
# A.U. chain MA: Image chain QB
# A.U. chain NA: Image chain RB
############################### Segments follow ################################
# Segment A begins
############################### Segment C begins ###############################
set C00 [atomselect $m1 "serial 1 to 1196"]
$C00 set segname C
set C00_data [ backup $C00 [ list chain x y z resid resname name ] ]
$C00 set chain C
$C00 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$C00 writepdb segtype_polymer_C_1_to_147.pdb
set C02 [atomselect $m1 "serial 1197 to 2817"]
$C02 set segname C
set C02_data [ backup $C02 [ list chain x y z resid resname name ] ]
$C02 set chain C
$C02 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$C02 writepdb segtype_polymer_C_157_to_367.pdb
set C04 [atomselect $m1 "serial 2818 to 3565"]
$C04 set segname C
set C04_data [ backup $C04 [ list chain x y z resid resname name ] ]
$C04 set chain C
$C04 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$C04 writepdb segtype_polymer_C_379_to_473.pdb
segment C {
    pdb segtype_polymer_C_1_to_147.pdb
    residue 148 GLU C
    residue 149 ASN C
    residue 150 GLN C
    residue 151 GLY C
    residue 152 ASN C
    residue 153 ARG C
    residue 154 SER C
    residue 155 ASN C
    residue 156 ASN C
    residue 156A GLY C
    pdb segtype_polymer_C_157_to_367.pdb
    residue 368 THR C
    residue 369 SER C
    residue 370 VAL C
    residue 371 GLN C
    residue 372 GLY C
    residue 373 SER C
    residue 374 ASN C
    residue 375 SER C
    residue 376 THR C
    residue 377 GLY C
    residue 378 SER C
    residue 378A GLY C
    pdb segtype_polymer_C_379_to_473.pdb
    mutate 469 ALA
}
################################ End segment C #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_C_1_to_147.pdb C
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_C_157_to_367.pdb C
# Subsegment [4] is a resolved run
coordpdb segtype_polymer_C_379_to_473.pdb C
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at A_GLU148 from A_ASN147
coord C 148 N {-268.67223 -168.66562 -46.16939}
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at A_THR368 from A_ASN367
coord C 368 N {-263.37298 -212.34323 -7.07846}
####################### Intra-segmental terminal patches #######################
patch CTER C:156
patch NTER C:157
delatom C 156A
patch CTER C:378
patch NTER C:379
delatom C 378A
############## Restoring A.U. state for all resolved subsegments ###############
restore $C00 [ list chain x y z resid resname name ]  $C00_data
restore $C02 [ list chain x y z resid resname name ]  $C02_data
restore $C04 [ list chain x y z resid resname name ]  $C04_data
################################ Segment C ends ################################
# Segment A ends
# Segment B begins
############################### Segment D begins ###############################
set D00 [atomselect $m1 "serial 3566 to 3769"]
$D00 set segname D
set D00_data [ backup $D00 [ list chain x y z resid resname name ] ]
$D00 set chain D
$D00 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$D00 writepdb segtype_polymer_D_7_to_36.pdb
set D02 [atomselect $m1 "serial 3770 to 4566"]
$D02 set segname D
set D02_data [ backup $D02 [ list chain x y z resid resname name ] ]
$D02 set chain D
$D02 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$D02 writepdb segtype_polymer_D_58_to_153.pdb
segment D {
    pdb segtype_polymer_D_7_to_36.pdb
    residue 37 ILE D
    residue 38 VAL D
    residue 39 GLN D
    residue 40 GLN D
    residue 41 GLN D
    residue 42 SER D
    residue 43 ASN D
    residue 44 LEU D
    residue 45 LEU D
    residue 46 ARG D
    residue 47 ALA D
    residue 48 PRO D
    residue 49 GLU D
    residue 50 ALA D
    residue 51 GLN D
    residue 52 GLN D
    residue 53 HSD D
    residue 54 LEU D
    residue 55 LEU D
    residue 56 LYS D
    residue 57 LEU D
    residue 57A GLY D
    pdb segtype_polymer_D_58_to_153.pdb
    mutate 94 THR
    mutate 48 ILE
}
################################ End segment D #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_D_7_to_36.pdb D
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_D_58_to_153.pdb D
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B_ILE37 from B_GLY36
coord D 37 N {-250.69074 -155.34257 18.97426}
####################### Intra-segmental terminal patches #######################
patch CTER D:57
patch NTER D:58
delatom D 57A
############## Restoring A.U. state for all resolved subsegments ###############
restore $D00 [ list chain x y z resid resname name ]  $D00_data
restore $D02 [ list chain x y z resid resname name ]  $D02_data
################################ Segment D ends ################################
# Segment B ends
# Segment AG01 begins
##################### Segment CG01 begins as image of AG01 #####################
set CG01 [atomselect $m1 "serial 11345 to 11427"]
$CG01 set segname CG01
set CG01_data [ backup $CG01 [ list chain x y z resid resname name ] ]
$CG01 set chain C
$CG01 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG01 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7]
$CG01 writepdb segtype_generic_CG01.pdb
segment CG01 {
    first none
    last none
    pdb segtype_generic_CG01.pdb
}
coordpdb segtype_generic_CG01.pdb CG01
######################## Restoring A.U. state for AG01 #########################
restore $CG01 [ list chain x y z resid resname name ]  $CG01_data
############################## Segment CG01 ends ###############################
# Segment AG01 ends
# Segment AG02 begins
##################### Segment CG02 begins as image of AG02 #####################
set CG02 [atomselect $m1 "serial 11428 to 11477"]
$CG02 set segname CG02
set CG02_data [ backup $CG02 [ list chain x y z resid resname name ] ]
$CG02 set chain C
$CG02 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG02 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4]
$CG02 writepdb segtype_generic_CG02.pdb
segment CG02 {
    first none
    last none
    pdb segtype_generic_CG02.pdb
}
coordpdb segtype_generic_CG02.pdb CG02
######################## Restoring A.U. state for AG02 #########################
restore $CG02 [ list chain x y z resid resname name ]  $CG02_data
############################## Segment CG02 ends ###############################
# Segment AG02 ends
# Segment AG03 begins
##################### Segment CG03 begins as image of AG03 #####################
set CG03 [atomselect $m1 "serial 11478 to 11505"]
$CG03 set segname CG03
set CG03_data [ backup $CG03 [ list chain x y z resid resname name ] ]
$CG03 set chain C
$CG03 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG03 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG03 writepdb segtype_generic_CG03.pdb
segment CG03 {
    first none
    last none
    pdb segtype_generic_CG03.pdb
}
coordpdb segtype_generic_CG03.pdb CG03
######################## Restoring A.U. state for AG03 #########################
restore $CG03 [ list chain x y z resid resname name ]  $CG03_data
############################## Segment CG03 ends ###############################
# Segment AG03 ends
# Segment AG04 begins
##################### Segment CG04 begins as image of AG04 #####################
set CG04 [atomselect $m1 "serial 11506 to 11577"]
$CG04 set segname CG04
set CG04_data [ backup $CG04 [ list chain x y z resid resname name ] ]
$CG04 set chain C
$CG04 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG04 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6]
$CG04 writepdb segtype_generic_CG04.pdb
segment CG04 {
    first none
    last none
    pdb segtype_generic_CG04.pdb
}
coordpdb segtype_generic_CG04.pdb CG04
######################## Restoring A.U. state for AG04 #########################
restore $CG04 [ list chain x y z resid resname name ]  $CG04_data
############################## Segment CG04 ends ###############################
# Segment AG04 ends
# Segment AG05 begins
##################### Segment CG05 begins as image of AG05 #####################
set CG05 [atomselect $m1 "serial 11578 to 11605"]
$CG05 set segname CG05
set CG05_data [ backup $CG05 [ list chain x y z resid resname name ] ]
$CG05 set chain C
$CG05 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG05 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG05 writepdb segtype_generic_CG05.pdb
segment CG05 {
    first none
    last none
    pdb segtype_generic_CG05.pdb
}
coordpdb segtype_generic_CG05.pdb CG05
######################## Restoring A.U. state for AG05 #########################
restore $CG05 [ list chain x y z resid resname name ]  $CG05_data
############################## Segment CG05 ends ###############################
# Segment AG05 ends
# Segment AG06 begins
##################### Segment CG06 begins as image of AG06 #####################
set CG06 [atomselect $m1 "serial 11606 to 11633"]
$CG06 set segname CG06
set CG06_data [ backup $CG06 [ list chain x y z resid resname name ] ]
$CG06 set chain C
$CG06 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG06 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG06 writepdb segtype_generic_CG06.pdb
segment CG06 {
    first none
    last none
    pdb segtype_generic_CG06.pdb
}
coordpdb segtype_generic_CG06.pdb CG06
######################## Restoring A.U. state for AG06 #########################
restore $CG06 [ list chain x y z resid resname name ]  $CG06_data
############################## Segment CG06 ends ###############################
# Segment AG06 ends
# Segment AG07 begins
##################### Segment CG07 begins as image of AG07 #####################
set CG07 [atomselect $m1 "serial 11634 to 11661"]
$CG07 set segname CG07
set CG07_data [ backup $CG07 [ list chain x y z resid resname name ] ]
$CG07 set chain C
$CG07 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG07 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG07 writepdb segtype_generic_CG07.pdb
segment CG07 {
    first none
    last none
    pdb segtype_generic_CG07.pdb
}
coordpdb segtype_generic_CG07.pdb CG07
######################## Restoring A.U. state for AG07 #########################
restore $CG07 [ list chain x y z resid resname name ]  $CG07_data
############################## Segment CG07 ends ###############################
# Segment AG07 ends
# Segment AG08 begins
##################### Segment CG08 begins as image of AG08 #####################
set CG08 [atomselect $m1 "serial 11662 to 11700"]
$CG08 set segname CG08
set CG08_data [ backup $CG08 [ list chain x y z resid resname name ] ]
$CG08 set chain C
$CG08 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG08 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3]
$CG08 writepdb segtype_generic_CG08.pdb
segment CG08 {
    first none
    last none
    pdb segtype_generic_CG08.pdb
}
coordpdb segtype_generic_CG08.pdb CG08
######################## Restoring A.U. state for AG08 #########################
restore $CG08 [ list chain x y z resid resname name ]  $CG08_data
############################## Segment CG08 ends ###############################
# Segment AG08 ends
# Segment AG09 begins
##################### Segment CG09 begins as image of AG09 #####################
set CG09 [atomselect $m1 "serial 11701 to 11761"]
$CG09 set segname CG09
set CG09_data [ backup $CG09 [ list chain x y z resid resname name ] ]
$CG09 set chain C
$CG09 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG09 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5]
$CG09 writepdb segtype_generic_CG09.pdb
segment CG09 {
    first none
    last none
    pdb segtype_generic_CG09.pdb
}
coordpdb segtype_generic_CG09.pdb CG09
######################## Restoring A.U. state for AG09 #########################
restore $CG09 [ list chain x y z resid resname name ]  $CG09_data
############################## Segment CG09 ends ###############################
# Segment AG09 ends
# Segment AG10 begins
##################### Segment CG10 begins as image of AG10 #####################
set CG10 [atomselect $m1 "serial 11762 to 11789"]
$CG10 set segname CG10
set CG10_data [ backup $CG10 [ list chain x y z resid resname name ] ]
$CG10 set chain C
$CG10 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG10 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG10 writepdb segtype_generic_CG10.pdb
segment CG10 {
    first none
    last none
    pdb segtype_generic_CG10.pdb
}
coordpdb segtype_generic_CG10.pdb CG10
######################## Restoring A.U. state for AG10 #########################
restore $CG10 [ list chain x y z resid resname name ]  $CG10_data
############################## Segment CG10 ends ###############################
# Segment AG10 ends
# Segment AG11 begins
##################### Segment CG11 begins as image of AG11 #####################
set CG11 [atomselect $m1 "serial 11790 to 11817"]
$CG11 set segname CG11
set CG11_data [ backup $CG11 [ list chain x y z resid resname name ] ]
$CG11 set chain C
$CG11 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG11 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG11 writepdb segtype_generic_CG11.pdb
segment CG11 {
    first none
    last none
    pdb segtype_generic_CG11.pdb
}
coordpdb segtype_generic_CG11.pdb CG11
######################## Restoring A.U. state for AG11 #########################
restore $CG11 [ list chain x y z resid resname name ]  $CG11_data
############################## Segment CG11 ends ###############################
# Segment AG11 ends
# Segment AG12 begins
##################### Segment CG12 begins as image of AG12 #####################
set CG12 [atomselect $m1 "serial 11818 to 11845"]
$CG12 set segname CG12
set CG12_data [ backup $CG12 [ list chain x y z resid resname name ] ]
$CG12 set chain C
$CG12 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG12 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG12 writepdb segtype_generic_CG12.pdb
segment CG12 {
    first none
    last none
    pdb segtype_generic_CG12.pdb
}
coordpdb segtype_generic_CG12.pdb CG12
######################## Restoring A.U. state for AG12 #########################
restore $CG12 [ list chain x y z resid resname name ]  $CG12_data
############################## Segment CG12 ends ###############################
# Segment AG12 ends
# Segment AG13 begins
##################### Segment CG13 begins as image of AG13 #####################
set CG13 [atomselect $m1 "serial 11846 to 11873"]
$CG13 set segname CG13
set CG13_data [ backup $CG13 [ list chain x y z resid resname name ] ]
$CG13 set chain C
$CG13 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$CG13 writepdb segtype_generic_CG13.pdb
segment CG13 {
    first none
    last none
    pdb segtype_generic_CG13.pdb
}
coordpdb segtype_generic_CG13.pdb CG13
######################## Restoring A.U. state for AG13 #########################
restore $CG13 [ list chain x y z resid resname name ]  $CG13_data
############################## Segment CG13 ends ###############################
# Segment AG13 ends
# Segment AG14 begins
##################### Segment CG14 begins as image of AG14 #####################
set CG14 [atomselect $m1 "serial 11874 to 11989"]
$CG14 set segname CG14
set CG14_data [ backup $CG14 [ list chain x y z resid resname name ] ]
$CG14 set chain C
$CG14 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG14 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10]
$CG14 writepdb segtype_generic_CG14.pdb
segment CG14 {
    first none
    last none
    pdb segtype_generic_CG14.pdb
}
coordpdb segtype_generic_CG14.pdb CG14
######################## Restoring A.U. state for AG14 #########################
restore $CG14 [ list chain x y z resid resname name ]  $CG14_data
############################## Segment CG14 ends ###############################
# Segment AG14 ends
# Segment AG15 begins
##################### Segment CG15 begins as image of AG15 #####################
set CG15 [atomselect $m1 "serial 11990 to 12003"]
$CG15 set segname CG15
set CG15_data [ backup $CG15 [ list chain x y z resid resname name ] ]
$CG15 set chain C
$CG15 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG15 set resid [list 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276]
$CG15 writepdb segtype_generic_CG15.pdb
segment CG15 {
    first none
    last none
    pdb segtype_generic_CG15.pdb
}
coordpdb segtype_generic_CG15.pdb CG15
######################## Restoring A.U. state for AG15 #########################
restore $CG15 [ list chain x y z resid resname name ]  $CG15_data
############################## Segment CG15 ends ###############################
# Segment AG15 ends
# Segment AG16 begins
##################### Segment CG16 begins as image of AG16 #####################
set CG16 [atomselect $m1 "serial 12004 to 12017"]
$CG16 set segname CG16
set CG16_data [ backup $CG16 [ list chain x y z resid resname name ] ]
$CG16 set chain C
$CG16 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG16 set resid [list 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839]
$CG16 writepdb segtype_generic_CG16.pdb
segment CG16 {
    first none
    last none
    pdb segtype_generic_CG16.pdb
}
coordpdb segtype_generic_CG16.pdb CG16
######################## Restoring A.U. state for AG16 #########################
restore $CG16 [ list chain x y z resid resname name ]  $CG16_data
############################## Segment CG16 ends ###############################
# Segment AG16 ends
# Segment AG17 begins
##################### Segment CG17 begins as image of AG17 #####################
set CG17 [atomselect $m1 "serial 12018 to 12031"]
$CG17 set segname CG17
set CG17_data [ backup $CG17 [ list chain x y z resid resname name ] ]
$CG17 set chain C
$CG17 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG17 set resid [list 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355]
$CG17 writepdb segtype_generic_CG17.pdb
segment CG17 {
    first none
    last none
    pdb segtype_generic_CG17.pdb
}
coordpdb segtype_generic_CG17.pdb CG17
######################## Restoring A.U. state for AG17 #########################
restore $CG17 [ list chain x y z resid resname name ]  $CG17_data
############################## Segment CG17 ends ###############################
# Segment AG17 ends
# Segment AG18 begins
##################### Segment CG18 begins as image of AG18 #####################
set CG18 [atomselect $m1 "serial 12032 to 12045"]
$CG18 set segname CG18
set CG18_data [ backup $CG18 [ list chain x y z resid resname name ] ]
$CG18 set chain C
$CG18 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$CG18 set resid [list 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133]
$CG18 writepdb segtype_generic_CG18.pdb
segment CG18 {
    first none
    last none
    pdb segtype_generic_CG18.pdb
}
coordpdb segtype_generic_CG18.pdb CG18
######################## Restoring A.U. state for AG18 #########################
restore $CG18 [ list chain x y z resid resname name ]  $CG18_data
############################## Segment CG18 ends ###############################
# Segment AG18 ends
# Segment BG01 begins
##################### Segment DG01 begins as image of BG01 #####################
set DG01 [atomselect $m1 "serial 12081 to 12094"]
$DG01 set segname DG01
set DG01_data [ backup $DG01 [ list chain x y z resid resname name ] ]
$DG01 set chain D
$DG01 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$DG01 set resid [list 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611]
$DG01 writepdb segtype_generic_DG01.pdb
segment DG01 {
    first none
    last none
    pdb segtype_generic_DG01.pdb
}
coordpdb segtype_generic_DG01.pdb DG01
######################## Restoring A.U. state for BG01 #########################
restore $DG01 [ list chain x y z resid resname name ]  $DG01_data
############################## Segment DG01 ends ###############################
# Segment BG01 ends
# Segment BG02 begins
##################### Segment DG02 begins as image of BG02 #####################
set DG02 [atomselect $m1 "serial 12095 to 12108"]
$DG02 set segname DG02
set DG02_data [ backup $DG02 [ list chain x y z resid resname name ] ]
$DG02 set chain D
$DG02 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$DG02 set resid [list 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618]
$DG02 writepdb segtype_generic_DG02.pdb
segment DG02 {
    first none
    last none
    pdb segtype_generic_DG02.pdb
}
coordpdb segtype_generic_DG02.pdb DG02
######################## Restoring A.U. state for BG02 #########################
restore $DG02 [ list chain x y z resid resname name ]  $DG02_data
############################## Segment DG02 ends ###############################
# Segment BG02 ends
# Segment BG03 begins
##################### Segment DG03 begins as image of BG03 #####################
set DG03 [atomselect $m1 "serial 12109 to 12122"]
$DG03 set segname DG03
set DG03_data [ backup $DG03 [ list chain x y z resid resname name ] ]
$DG03 set chain D
$DG03 move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$DG03 set resid [list 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637]
$DG03 writepdb segtype_generic_DG03.pdb
segment DG03 {
    first none
    last none
    pdb segtype_generic_DG03.pdb
}
coordpdb segtype_generic_DG03.pdb DG03
######################## Restoring A.U. state for BG03 #########################
restore $DG03 [ list chain x y z resid resname name ]  $DG03_data
############################## Segment DG03 ends ###############################
# Segment BG03 ends
# Segment Y begins
####################### Segment HB begins as image of Y ########################
set HB [atomselect $m1 "serial 12046 to 12050"]
$HB set segname HB
set HB_data [ backup $HB [ list chain x y z resid resname name ] ]
$HB set chain HB
$HB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$HB set resid [list 601 601 601 601 601]
$HB writepdb segtype_generic_HB.pdb
segment HB {
    first none
    last none
    pdb segtype_generic_HB.pdb
}
coordpdb segtype_generic_HB.pdb HB
########################## Restoring A.U. state for Y ##########################
restore $HB [ list chain x y z resid resname name ]  $HB_data
############################### Segment HB ends ################################
# Segment Y ends
# Segment Z begins
####################### Segment IB begins as image of Z ########################
set IB [atomselect $m1 "serial 12051 to 12055"]
$IB set segname IB
set IB_data [ backup $IB [ list chain x y z resid resname name ] ]
$IB set chain IB
$IB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$IB set resid [list 602 602 602 602 602]
$IB writepdb segtype_generic_IB.pdb
segment IB {
    first none
    last none
    pdb segtype_generic_IB.pdb
}
coordpdb segtype_generic_IB.pdb IB
########################## Restoring A.U. state for Z ##########################
restore $IB [ list chain x y z resid resname name ]  $IB_data
############################### Segment IB ends ################################
# Segment Z ends
# Segment AA begins
####################### Segment JB begins as image of AA #######################
set JB [atomselect $m1 "serial 12056 to 12060"]
$JB set segname JB
set JB_data [ backup $JB [ list chain x y z resid resname name ] ]
$JB set chain JB
$JB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$JB set resid [list 603 603 603 603 603]
$JB writepdb segtype_generic_JB.pdb
segment JB {
    first none
    last none
    pdb segtype_generic_JB.pdb
}
coordpdb segtype_generic_JB.pdb JB
######################### Restoring A.U. state for AA ##########################
restore $JB [ list chain x y z resid resname name ]  $JB_data
############################### Segment JB ends ################################
# Segment AA ends
# Segment BA begins
####################### Segment KB begins as image of BA #######################
set KB [atomselect $m1 "serial 12061 to 12065"]
$KB set segname KB
set KB_data [ backup $KB [ list chain x y z resid resname name ] ]
$KB set chain KB
$KB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$KB set resid [list 604 604 604 604 604]
$KB writepdb segtype_generic_KB.pdb
segment KB {
    first none
    last none
    pdb segtype_generic_KB.pdb
}
coordpdb segtype_generic_KB.pdb KB
######################### Restoring A.U. state for BA ##########################
restore $KB [ list chain x y z resid resname name ]  $KB_data
############################### Segment KB ends ################################
# Segment BA ends
# Segment CA begins
####################### Segment LB begins as image of CA #######################
set LB [atomselect $m1 "serial 12066 to 12070"]
$LB set segname LB
set LB_data [ backup $LB [ list chain x y z resid resname name ] ]
$LB set chain LB
$LB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$LB set resid [list 605 605 605 605 605]
$LB writepdb segtype_generic_LB.pdb
segment LB {
    first none
    last none
    pdb segtype_generic_LB.pdb
}
coordpdb segtype_generic_LB.pdb LB
######################### Restoring A.U. state for CA ##########################
restore $LB [ list chain x y z resid resname name ]  $LB_data
############################### Segment LB ends ################################
# Segment CA ends
# Segment DA begins
####################### Segment MB begins as image of DA #######################
set MB [atomselect $m1 "serial 12071 to 12075"]
$MB set segname MB
set MB_data [ backup $MB [ list chain x y z resid resname name ] ]
$MB set chain MB
$MB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$MB set resid [list 606 606 606 606 606]
$MB writepdb segtype_generic_MB.pdb
segment MB {
    first none
    last none
    pdb segtype_generic_MB.pdb
}
coordpdb segtype_generic_MB.pdb MB
######################### Restoring A.U. state for DA ##########################
restore $MB [ list chain x y z resid resname name ]  $MB_data
############################### Segment MB ends ################################
# Segment DA ends
# Segment EA begins
####################### Segment NB begins as image of EA #######################
set NB [atomselect $m1 "serial 12076 to 12080"]
$NB set segname NB
set NB_data [ backup $NB [ list chain x y z resid resname name ] ]
$NB set chain NB
$NB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$NB set resid [list 607 607 607 607 607]
$NB writepdb segtype_generic_NB.pdb
segment NB {
    first none
    last none
    pdb segtype_generic_NB.pdb
}
coordpdb segtype_generic_NB.pdb NB
######################### Restoring A.U. state for EA ##########################
restore $NB [ list chain x y z resid resname name ]  $NB_data
############################### Segment NB ends ################################
# Segment EA ends
# Segment IA begins
####################### Segment OB begins as image of IA #######################
set OB [atomselect $m1 "serial 12123 to 12127"]
$OB set segname OB
set OB_data [ backup $OB [ list chain x y z resid resname name ] ]
$OB set chain OB
$OB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$OB set resid [list 701 701 701 701 701]
$OB writepdb segtype_generic_OB.pdb
segment OB {
    first none
    last none
    pdb segtype_generic_OB.pdb
}
coordpdb segtype_generic_OB.pdb OB
######################### Restoring A.U. state for IA ##########################
restore $OB [ list chain x y z resid resname name ]  $OB_data
############################### Segment OB ends ################################
# Segment IA ends
# Segment JA begins
####################### Segment PB begins as image of JA #######################
set PB [atomselect $m1 "serial 12128 to 12132"]
$PB set segname PB
set PB_data [ backup $PB [ list chain x y z resid resname name ] ]
$PB set chain PB
$PB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$PB set resid [list 702 702 702 702 702]
$PB writepdb segtype_generic_PB.pdb
segment PB {
    first none
    last none
    pdb segtype_generic_PB.pdb
}
coordpdb segtype_generic_PB.pdb PB
######################### Restoring A.U. state for JA ##########################
restore $PB [ list chain x y z resid resname name ]  $PB_data
############################### Segment PB ends ################################
# Segment JA ends
# Segment MA begins
####################### Segment QB begins as image of MA #######################
set QB [atomselect $m1 "serial 12152 to 12168"]
$QB set segname QB
set QB_data [ backup $QB [ list chain x y z resid resname name ] ]
$QB set chain QB
$QB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$QB set resid [list 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717]
$QB writepdb segtype_generic_QB.pdb
segment QB {
    first none
    last none
    pdb segtype_generic_QB.pdb
}
coordpdb segtype_generic_QB.pdb QB
######################### Restoring A.U. state for MA ##########################
restore $QB [ list chain x y z resid resname name ]  $QB_data
############################### Segment QB ends ################################
# Segment MA ends
# Segment NA begins
####################### Segment RB begins as image of NA #######################
set RB [atomselect $m1 "serial 12169 12170"]
$RB set segname RB
set RB_data [ backup $RB [ list chain x y z resid resname name ] ]
$RB set chain RB
$RB move { { -0.5 -0.8660254038 0.0 -515.56  } { 0.8660254038 -0.5 0.0 0.0  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$RB set resid [list 801 802]
$RB writepdb segtype_generic_RB.pdb
segment RB {
    first none
    last none
    pdb segtype_generic_RB.pdb
}
coordpdb segtype_generic_RB.pdb RB
######################### Restoring A.U. state for NA ##########################
restore $RB [ list chain x y z resid resname name ]  $RB_data
############################### Segment RB ends ################################
# Segment NA ends
############################# DISU patches follow ##############################
patch DISU C:24 C:44
patch DISU C:89 C:175
patch DISU C:96 C:166
patch DISU C:101 C:119
patch DISU C:188 C:217
patch DISU C:198 C:209
patch DISU C:266 C:300
patch DISU C:347 C:413
patch DISU C:354 C:386
patch DISU D:87 D:93
############################# LINK patches follow ##############################
patch NGLA C:58 CG01:1
patch NGLB C:103 CG18:1133
patch NGLA C:107 CG02:1
patch NGLA C:118 CG09:1
patch NGLA C:122 CG10:1
patch NGLA C:167 CG11:1
patch NGLA C:204 CG03:1
patch NGLA C:232 CG04:1
patch NGLA C:246 CG15:1276
patch NGLA C:265 CG12:1
patch NGLB C:271 CG13:1
patch NGLA C:301 CG14:1
patch NGLB C:308 CG16:1839
patch NGLA C:324 CG17:1355
patch NGLA C:332 CG05:1
patch NGLA C:355 CG06:1
patch NGLB C:361 CG07:1
patch NGLA C:416 CG08:1
patch NGLB D:100 DG01:1611
patch NGLA D:107 DG02:1618
patch NGLA D:126 DG03:1637
patch 14ab CG01:1 CG01:2
patch 14ab CG01:2 CG01:3
patch 16BT CG01:3 CG01:4
patch 13aa CG01:3 CG01:7
patch 13ba CG01:4 CG01:5
patch 16AT CG01:4 CG01:6
patch 14ab CG02:1 CG02:2
patch 14ab CG02:2 CG02:3
patch 13ba CG02:3 CG02:4
patch 14ab CG03:1 CG03:2
patch 14bb CG04:1 CG04:2
patch 14aa CG04:2 CG04:3
patch 13ba CG04:3 CG04:4
patch 16AT CG04:3 CG04:6
patch 12ba CG04:4 CG04:5
patch 14bb CG05:1 CG05:2
patch 14bb CG06:1 CG06:2
patch 14bb CG07:1 CG07:2
patch 14ab CG08:1 CG08:2
patch 14bb CG08:2 CG08:3
patch 14ab CG09:1 CG09:2
patch 14ab CG09:2 CG09:3
patch 13ba CG09:3 CG09:4
patch 16AT CG09:3 CG09:5
patch 14ab CG10:1 CG10:2
patch 14ab CG11:1 CG11:2
patch 14ab CG12:1 CG12:2
patch 14bb CG13:1 CG13:2
patch 14bb CG14:1 CG14:2
patch 14bb CG14:2 CG14:3
patch 13aa CG14:3 CG14:4
patch 16AT CG14:3 CG14:7
patch 12aa CG14:4 CG14:5
patch 12ab CG14:5 CG14:6
patch 13bb CG14:7 CG14:8
patch 16AT CG14:7 CG14:10
patch 12bb CG14:8 CG14:9
############################### Transform 1 ends ###############################
############################## Transform 2 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
# A.U. chain A: Image chain E
# A.U. chain B: Image chain F
# A.U. chain AG01: Image chain EG01
# A.U. chain AG02: Image chain EG02
# A.U. chain AG03: Image chain EG03
# A.U. chain AG04: Image chain EG04
# A.U. chain AG05: Image chain EG05
# A.U. chain AG06: Image chain EG06
# A.U. chain AG07: Image chain EG07
# A.U. chain AG08: Image chain EG08
# A.U. chain AG09: Image chain EG09
# A.U. chain AG10: Image chain EG10
# A.U. chain AG11: Image chain EG11
# A.U. chain AG12: Image chain EG12
# A.U. chain AG13: Image chain EG13
# A.U. chain AG14: Image chain EG14
# A.U. chain AG15: Image chain EG15
# A.U. chain AG16: Image chain EG16
# A.U. chain AG17: Image chain EG17
# A.U. chain AG18: Image chain EG18
# A.U. chain BG01: Image chain FG01
# A.U. chain BG02: Image chain FG02
# A.U. chain BG03: Image chain FG03
# A.U. chain Y: Image chain NC
# A.U. chain Z: Image chain OC
# A.U. chain AA: Image chain PC
# A.U. chain BA: Image chain QC
# A.U. chain CA: Image chain RC
# A.U. chain DA: Image chain SC
# A.U. chain EA: Image chain TC
# A.U. chain IA: Image chain UC
# A.U. chain JA: Image chain VC
# A.U. chain MA: Image chain WC
# A.U. chain NA: Image chain XC
############################### Segments follow ################################
# Segment A begins
############################### Segment E begins ###############################
set E00 [atomselect $m1 "serial 1 to 1196"]
$E00 set segname E
set E00_data [ backup $E00 [ list chain x y z resid resname name ] ]
$E00 set chain E
$E00 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$E00 writepdb segtype_polymer_E_1_to_147.pdb
set E02 [atomselect $m1 "serial 1197 to 2817"]
$E02 set segname E
set E02_data [ backup $E02 [ list chain x y z resid resname name ] ]
$E02 set chain E
$E02 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$E02 writepdb segtype_polymer_E_157_to_367.pdb
set E04 [atomselect $m1 "serial 2818 to 3565"]
$E04 set segname E
set E04_data [ backup $E04 [ list chain x y z resid resname name ] ]
$E04 set chain E
$E04 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$E04 writepdb segtype_polymer_E_379_to_473.pdb
segment E {
    pdb segtype_polymer_E_1_to_147.pdb
    residue 148 GLU E
    residue 149 ASN E
    residue 150 GLN E
    residue 151 GLY E
    residue 152 ASN E
    residue 153 ARG E
    residue 154 SER E
    residue 155 ASN E
    residue 156 ASN E
    residue 156A GLY E
    pdb segtype_polymer_E_157_to_367.pdb
    residue 368 THR E
    residue 369 SER E
    residue 370 VAL E
    residue 371 GLN E
    residue 372 GLY E
    residue 373 SER E
    residue 374 ASN E
    residue 375 SER E
    residue 376 THR E
    residue 377 GLY E
    residue 378 SER E
    residue 378A GLY E
    pdb segtype_polymer_E_379_to_473.pdb
    mutate 469 ALA
}
################################ End segment E #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_E_1_to_147.pdb E
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_E_157_to_367.pdb E
# Subsegment [4] is a resolved run
coordpdb segtype_polymer_E_379_to_473.pdb E
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at A_GLU148 from A_ASN147
coord E 148 N {-235.15518 -148.34417 -46.16939}
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at A_THR368 from A_ASN367
coord E 368 N {-199.97888 -121.91607 -7.07846}
####################### Intra-segmental terminal patches #######################
patch CTER E:156
patch NTER E:157
delatom E 156A
patch CTER E:378
patch NTER E:379
delatom E 378A
############## Restoring A.U. state for all resolved subsegments ###############
restore $E00 [ list chain x y z resid resname name ]  $E00_data
restore $E02 [ list chain x y z resid resname name ]  $E02_data
restore $E04 [ list chain x y z resid resname name ]  $E04_data
################################ Segment E ends ################################
# Segment A ends
# Segment B begins
############################### Segment F begins ###############################
set F00 [atomselect $m1 "serial 3566 to 3769"]
$F00 set segname F
set F00_data [ backup $F00 [ list chain x y z resid resname name ] ]
$F00 set chain F
$F00 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$F00 writepdb segtype_polymer_F_7_to_36.pdb
set F02 [atomselect $m1 "serial 3770 to 4566"]
$F02 set segname F
set F02_data [ backup $F02 [ list chain x y z resid resname name ] ]
$F02 set chain F
$F02 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$F02 writepdb segtype_polymer_F_58_to_153.pdb
segment F {
    pdb segtype_polymer_F_7_to_36.pdb
    residue 37 ILE F
    residue 38 VAL F
    residue 39 GLN F
    residue 40 GLN F
    residue 41 GLN F
    residue 42 SER F
    residue 43 ASN F
    residue 44 LEU F
    residue 45 LEU F
    residue 46 ARG F
    residue 47 ALA F
    residue 48 PRO F
    residue 49 GLU F
    residue 50 ALA F
    residue 51 GLN F
    residue 52 GLN F
    residue 53 HSD F
    residue 54 LEU F
    residue 55 LEU F
    residue 56 LYS F
    residue 57 LEU F
    residue 57A GLY F
    pdb segtype_polymer_F_58_to_153.pdb
    mutate 94 THR
    mutate 48 ILE
}
################################ End segment F #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_F_7_to_36.pdb F
# Subsegment [2] is a resolved run
coordpdb segtype_polymer_F_58_to_153.pdb F
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B_ILE37 from B_GLY36
coord F 37 N {-255.68402 -139.43326 18.97426}
####################### Intra-segmental terminal patches #######################
patch CTER F:57
patch NTER F:58
delatom F 57A
############## Restoring A.U. state for all resolved subsegments ###############
restore $F00 [ list chain x y z resid resname name ]  $F00_data
restore $F02 [ list chain x y z resid resname name ]  $F02_data
################################ Segment F ends ################################
# Segment B ends
# Segment AG01 begins
##################### Segment EG01 begins as image of AG01 #####################
set EG01 [atomselect $m1 "serial 11345 to 11427"]
$EG01 set segname EG01
set EG01_data [ backup $EG01 [ list chain x y z resid resname name ] ]
$EG01 set chain E
$EG01 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG01 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7]
$EG01 writepdb segtype_generic_EG01.pdb
segment EG01 {
    first none
    last none
    pdb segtype_generic_EG01.pdb
}
coordpdb segtype_generic_EG01.pdb EG01
######################## Restoring A.U. state for AG01 #########################
restore $EG01 [ list chain x y z resid resname name ]  $EG01_data
############################## Segment EG01 ends ###############################
# Segment AG01 ends
# Segment AG02 begins
##################### Segment EG02 begins as image of AG02 #####################
set EG02 [atomselect $m1 "serial 11428 to 11477"]
$EG02 set segname EG02
set EG02_data [ backup $EG02 [ list chain x y z resid resname name ] ]
$EG02 set chain E
$EG02 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG02 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4]
$EG02 writepdb segtype_generic_EG02.pdb
segment EG02 {
    first none
    last none
    pdb segtype_generic_EG02.pdb
}
coordpdb segtype_generic_EG02.pdb EG02
######################## Restoring A.U. state for AG02 #########################
restore $EG02 [ list chain x y z resid resname name ]  $EG02_data
############################## Segment EG02 ends ###############################
# Segment AG02 ends
# Segment AG03 begins
##################### Segment EG03 begins as image of AG03 #####################
set EG03 [atomselect $m1 "serial 11478 to 11505"]
$EG03 set segname EG03
set EG03_data [ backup $EG03 [ list chain x y z resid resname name ] ]
$EG03 set chain E
$EG03 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG03 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG03 writepdb segtype_generic_EG03.pdb
segment EG03 {
    first none
    last none
    pdb segtype_generic_EG03.pdb
}
coordpdb segtype_generic_EG03.pdb EG03
######################## Restoring A.U. state for AG03 #########################
restore $EG03 [ list chain x y z resid resname name ]  $EG03_data
############################## Segment EG03 ends ###############################
# Segment AG03 ends
# Segment AG04 begins
##################### Segment EG04 begins as image of AG04 #####################
set EG04 [atomselect $m1 "serial 11506 to 11577"]
$EG04 set segname EG04
set EG04_data [ backup $EG04 [ list chain x y z resid resname name ] ]
$EG04 set chain E
$EG04 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG04 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6]
$EG04 writepdb segtype_generic_EG04.pdb
segment EG04 {
    first none
    last none
    pdb segtype_generic_EG04.pdb
}
coordpdb segtype_generic_EG04.pdb EG04
######################## Restoring A.U. state for AG04 #########################
restore $EG04 [ list chain x y z resid resname name ]  $EG04_data
############################## Segment EG04 ends ###############################
# Segment AG04 ends
# Segment AG05 begins
##################### Segment EG05 begins as image of AG05 #####################
set EG05 [atomselect $m1 "serial 11578 to 11605"]
$EG05 set segname EG05
set EG05_data [ backup $EG05 [ list chain x y z resid resname name ] ]
$EG05 set chain E
$EG05 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG05 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG05 writepdb segtype_generic_EG05.pdb
segment EG05 {
    first none
    last none
    pdb segtype_generic_EG05.pdb
}
coordpdb segtype_generic_EG05.pdb EG05
######################## Restoring A.U. state for AG05 #########################
restore $EG05 [ list chain x y z resid resname name ]  $EG05_data
############################## Segment EG05 ends ###############################
# Segment AG05 ends
# Segment AG06 begins
##################### Segment EG06 begins as image of AG06 #####################
set EG06 [atomselect $m1 "serial 11606 to 11633"]
$EG06 set segname EG06
set EG06_data [ backup $EG06 [ list chain x y z resid resname name ] ]
$EG06 set chain E
$EG06 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG06 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG06 writepdb segtype_generic_EG06.pdb
segment EG06 {
    first none
    last none
    pdb segtype_generic_EG06.pdb
}
coordpdb segtype_generic_EG06.pdb EG06
######################## Restoring A.U. state for AG06 #########################
restore $EG06 [ list chain x y z resid resname name ]  $EG06_data
############################## Segment EG06 ends ###############################
# Segment AG06 ends
# Segment AG07 begins
##################### Segment EG07 begins as image of AG07 #####################
set EG07 [atomselect $m1 "serial 11634 to 11661"]
$EG07 set segname EG07
set EG07_data [ backup $EG07 [ list chain x y z resid resname name ] ]
$EG07 set chain E
$EG07 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG07 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG07 writepdb segtype_generic_EG07.pdb
segment EG07 {
    first none
    last none
    pdb segtype_generic_EG07.pdb
}
coordpdb segtype_generic_EG07.pdb EG07
######################## Restoring A.U. state for AG07 #########################
restore $EG07 [ list chain x y z resid resname name ]  $EG07_data
############################## Segment EG07 ends ###############################
# Segment AG07 ends
# Segment AG08 begins
##################### Segment EG08 begins as image of AG08 #####################
set EG08 [atomselect $m1 "serial 11662 to 11700"]
$EG08 set segname EG08
set EG08_data [ backup $EG08 [ list chain x y z resid resname name ] ]
$EG08 set chain E
$EG08 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG08 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3]
$EG08 writepdb segtype_generic_EG08.pdb
segment EG08 {
    first none
    last none
    pdb segtype_generic_EG08.pdb
}
coordpdb segtype_generic_EG08.pdb EG08
######################## Restoring A.U. state for AG08 #########################
restore $EG08 [ list chain x y z resid resname name ]  $EG08_data
############################## Segment EG08 ends ###############################
# Segment AG08 ends
# Segment AG09 begins
##################### Segment EG09 begins as image of AG09 #####################
set EG09 [atomselect $m1 "serial 11701 to 11761"]
$EG09 set segname EG09
set EG09_data [ backup $EG09 [ list chain x y z resid resname name ] ]
$EG09 set chain E
$EG09 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG09 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5]
$EG09 writepdb segtype_generic_EG09.pdb
segment EG09 {
    first none
    last none
    pdb segtype_generic_EG09.pdb
}
coordpdb segtype_generic_EG09.pdb EG09
######################## Restoring A.U. state for AG09 #########################
restore $EG09 [ list chain x y z resid resname name ]  $EG09_data
############################## Segment EG09 ends ###############################
# Segment AG09 ends
# Segment AG10 begins
##################### Segment EG10 begins as image of AG10 #####################
set EG10 [atomselect $m1 "serial 11762 to 11789"]
$EG10 set segname EG10
set EG10_data [ backup $EG10 [ list chain x y z resid resname name ] ]
$EG10 set chain E
$EG10 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG10 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG10 writepdb segtype_generic_EG10.pdb
segment EG10 {
    first none
    last none
    pdb segtype_generic_EG10.pdb
}
coordpdb segtype_generic_EG10.pdb EG10
######################## Restoring A.U. state for AG10 #########################
restore $EG10 [ list chain x y z resid resname name ]  $EG10_data
############################## Segment EG10 ends ###############################
# Segment AG10 ends
# Segment AG11 begins
##################### Segment EG11 begins as image of AG11 #####################
set EG11 [atomselect $m1 "serial 11790 to 11817"]
$EG11 set segname EG11
set EG11_data [ backup $EG11 [ list chain x y z resid resname name ] ]
$EG11 set chain E
$EG11 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG11 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG11 writepdb segtype_generic_EG11.pdb
segment EG11 {
    first none
    last none
    pdb segtype_generic_EG11.pdb
}
coordpdb segtype_generic_EG11.pdb EG11
######################## Restoring A.U. state for AG11 #########################
restore $EG11 [ list chain x y z resid resname name ]  $EG11_data
############################## Segment EG11 ends ###############################
# Segment AG11 ends
# Segment AG12 begins
##################### Segment EG12 begins as image of AG12 #####################
set EG12 [atomselect $m1 "serial 11818 to 11845"]
$EG12 set segname EG12
set EG12_data [ backup $EG12 [ list chain x y z resid resname name ] ]
$EG12 set chain E
$EG12 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG12 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG12 writepdb segtype_generic_EG12.pdb
segment EG12 {
    first none
    last none
    pdb segtype_generic_EG12.pdb
}
coordpdb segtype_generic_EG12.pdb EG12
######################## Restoring A.U. state for AG12 #########################
restore $EG12 [ list chain x y z resid resname name ]  $EG12_data
############################## Segment EG12 ends ###############################
# Segment AG12 ends
# Segment AG13 begins
##################### Segment EG13 begins as image of AG13 #####################
set EG13 [atomselect $m1 "serial 11846 to 11873"]
$EG13 set segname EG13
set EG13_data [ backup $EG13 [ list chain x y z resid resname name ] ]
$EG13 set chain E
$EG13 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$EG13 writepdb segtype_generic_EG13.pdb
segment EG13 {
    first none
    last none
    pdb segtype_generic_EG13.pdb
}
coordpdb segtype_generic_EG13.pdb EG13
######################## Restoring A.U. state for AG13 #########################
restore $EG13 [ list chain x y z resid resname name ]  $EG13_data
############################## Segment EG13 ends ###############################
# Segment AG13 ends
# Segment AG14 begins
##################### Segment EG14 begins as image of AG14 #####################
set EG14 [atomselect $m1 "serial 11874 to 11989"]
$EG14 set segname EG14
set EG14_data [ backup $EG14 [ list chain x y z resid resname name ] ]
$EG14 set chain E
$EG14 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG14 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10]
$EG14 writepdb segtype_generic_EG14.pdb
segment EG14 {
    first none
    last none
    pdb segtype_generic_EG14.pdb
}
coordpdb segtype_generic_EG14.pdb EG14
######################## Restoring A.U. state for AG14 #########################
restore $EG14 [ list chain x y z resid resname name ]  $EG14_data
############################## Segment EG14 ends ###############################
# Segment AG14 ends
# Segment AG15 begins
##################### Segment EG15 begins as image of AG15 #####################
set EG15 [atomselect $m1 "serial 11990 to 12003"]
$EG15 set segname EG15
set EG15_data [ backup $EG15 [ list chain x y z resid resname name ] ]
$EG15 set chain E
$EG15 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG15 set resid [list 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276 1276]
$EG15 writepdb segtype_generic_EG15.pdb
segment EG15 {
    first none
    last none
    pdb segtype_generic_EG15.pdb
}
coordpdb segtype_generic_EG15.pdb EG15
######################## Restoring A.U. state for AG15 #########################
restore $EG15 [ list chain x y z resid resname name ]  $EG15_data
############################## Segment EG15 ends ###############################
# Segment AG15 ends
# Segment AG16 begins
##################### Segment EG16 begins as image of AG16 #####################
set EG16 [atomselect $m1 "serial 12004 to 12017"]
$EG16 set segname EG16
set EG16_data [ backup $EG16 [ list chain x y z resid resname name ] ]
$EG16 set chain E
$EG16 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG16 set resid [list 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839 1839]
$EG16 writepdb segtype_generic_EG16.pdb
segment EG16 {
    first none
    last none
    pdb segtype_generic_EG16.pdb
}
coordpdb segtype_generic_EG16.pdb EG16
######################## Restoring A.U. state for AG16 #########################
restore $EG16 [ list chain x y z resid resname name ]  $EG16_data
############################## Segment EG16 ends ###############################
# Segment AG16 ends
# Segment AG17 begins
##################### Segment EG17 begins as image of AG17 #####################
set EG17 [atomselect $m1 "serial 12018 to 12031"]
$EG17 set segname EG17
set EG17_data [ backup $EG17 [ list chain x y z resid resname name ] ]
$EG17 set chain E
$EG17 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG17 set resid [list 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355 1355]
$EG17 writepdb segtype_generic_EG17.pdb
segment EG17 {
    first none
    last none
    pdb segtype_generic_EG17.pdb
}
coordpdb segtype_generic_EG17.pdb EG17
######################## Restoring A.U. state for AG17 #########################
restore $EG17 [ list chain x y z resid resname name ]  $EG17_data
############################## Segment EG17 ends ###############################
# Segment AG17 ends
# Segment AG18 begins
##################### Segment EG18 begins as image of AG18 #####################
set EG18 [atomselect $m1 "serial 12032 to 12045"]
$EG18 set segname EG18
set EG18_data [ backup $EG18 [ list chain x y z resid resname name ] ]
$EG18 set chain E
$EG18 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$EG18 set resid [list 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133 1133]
$EG18 writepdb segtype_generic_EG18.pdb
segment EG18 {
    first none
    last none
    pdb segtype_generic_EG18.pdb
}
coordpdb segtype_generic_EG18.pdb EG18
######################## Restoring A.U. state for AG18 #########################
restore $EG18 [ list chain x y z resid resname name ]  $EG18_data
############################## Segment EG18 ends ###############################
# Segment AG18 ends
# Segment BG01 begins
##################### Segment FG01 begins as image of BG01 #####################
set FG01 [atomselect $m1 "serial 12081 to 12094"]
$FG01 set segname FG01
set FG01_data [ backup $FG01 [ list chain x y z resid resname name ] ]
$FG01 set chain F
$FG01 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$FG01 set resid [list 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611 1611]
$FG01 writepdb segtype_generic_FG01.pdb
segment FG01 {
    first none
    last none
    pdb segtype_generic_FG01.pdb
}
coordpdb segtype_generic_FG01.pdb FG01
######################## Restoring A.U. state for BG01 #########################
restore $FG01 [ list chain x y z resid resname name ]  $FG01_data
############################## Segment FG01 ends ###############################
# Segment BG01 ends
# Segment BG02 begins
##################### Segment FG02 begins as image of BG02 #####################
set FG02 [atomselect $m1 "serial 12095 to 12108"]
$FG02 set segname FG02
set FG02_data [ backup $FG02 [ list chain x y z resid resname name ] ]
$FG02 set chain F
$FG02 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$FG02 set resid [list 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618 1618]
$FG02 writepdb segtype_generic_FG02.pdb
segment FG02 {
    first none
    last none
    pdb segtype_generic_FG02.pdb
}
coordpdb segtype_generic_FG02.pdb FG02
######################## Restoring A.U. state for BG02 #########################
restore $FG02 [ list chain x y z resid resname name ]  $FG02_data
############################## Segment FG02 ends ###############################
# Segment BG02 ends
# Segment BG03 begins
##################### Segment FG03 begins as image of BG03 #####################
set FG03 [atomselect $m1 "serial 12109 to 12122"]
$FG03 set segname FG03
set FG03_data [ backup $FG03 [ list chain x y z resid resname name ] ]
$FG03 set chain F
$FG03 move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$FG03 set resid [list 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637 1637]
$FG03 writepdb segtype_generic_FG03.pdb
segment FG03 {
    first none
    last none
    pdb segtype_generic_FG03.pdb
}
coordpdb segtype_generic_FG03.pdb FG03
######################## Restoring A.U. state for BG03 #########################
restore $FG03 [ list chain x y z resid resname name ]  $FG03_data
############################## Segment FG03 ends ###############################
# Segment BG03 ends
# Segment Y begins
####################### Segment NC begins as image of Y ########################
set NC [atomselect $m1 "serial 12046 to 12050"]
$NC set segname NC
set NC_data [ backup $NC [ list chain x y z resid resname name ] ]
$NC set chain NC
$NC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$NC set resid [list 601 601 601 601 601]
$NC writepdb segtype_generic_NC.pdb
segment NC {
    first none
    last none
    pdb segtype_generic_NC.pdb
}
coordpdb segtype_generic_NC.pdb NC
########################## Restoring A.U. state for Y ##########################
restore $NC [ list chain x y z resid resname name ]  $NC_data
############################### Segment NC ends ################################
# Segment Y ends
# Segment Z begins
####################### Segment OC begins as image of Z ########################
set OC [atomselect $m1 "serial 12051 to 12055"]
$OC set segname OC
set OC_data [ backup $OC [ list chain x y z resid resname name ] ]
$OC set chain OC
$OC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$OC set resid [list 602 602 602 602 602]
$OC writepdb segtype_generic_OC.pdb
segment OC {
    first none
    last none
    pdb segtype_generic_OC.pdb
}
coordpdb segtype_generic_OC.pdb OC
########################## Restoring A.U. state for Z ##########################
restore $OC [ list chain x y z resid resname name ]  $OC_data
############################### Segment OC ends ################################
# Segment Z ends
# Segment AA begins
####################### Segment PC begins as image of AA #######################
set PC [atomselect $m1 "serial 12056 to 12060"]
$PC set segname PC
set PC_data [ backup $PC [ list chain x y z resid resname name ] ]
$PC set chain PC
$PC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$PC set resid [list 603 603 603 603 603]
$PC writepdb segtype_generic_PC.pdb
segment PC {
    first none
    last none
    pdb segtype_generic_PC.pdb
}
coordpdb segtype_generic_PC.pdb PC
######################### Restoring A.U. state for AA ##########################
restore $PC [ list chain x y z resid resname name ]  $PC_data
############################### Segment PC ends ################################
# Segment AA ends
# Segment BA begins
####################### Segment QC begins as image of BA #######################
set QC [atomselect $m1 "serial 12061 to 12065"]
$QC set segname QC
set QC_data [ backup $QC [ list chain x y z resid resname name ] ]
$QC set chain QC
$QC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$QC set resid [list 604 604 604 604 604]
$QC writepdb segtype_generic_QC.pdb
segment QC {
    first none
    last none
    pdb segtype_generic_QC.pdb
}
coordpdb segtype_generic_QC.pdb QC
######################### Restoring A.U. state for BA ##########################
restore $QC [ list chain x y z resid resname name ]  $QC_data
############################### Segment QC ends ################################
# Segment BA ends
# Segment CA begins
####################### Segment RC begins as image of CA #######################
set RC [atomselect $m1 "serial 12066 to 12070"]
$RC set segname RC
set RC_data [ backup $RC [ list chain x y z resid resname name ] ]
$RC set chain RC
$RC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$RC set resid [list 605 605 605 605 605]
$RC writepdb segtype_generic_RC.pdb
segment RC {
    first none
    last none
    pdb segtype_generic_RC.pdb
}
coordpdb segtype_generic_RC.pdb RC
######################### Restoring A.U. state for CA ##########################
restore $RC [ list chain x y z resid resname name ]  $RC_data
############################### Segment RC ends ################################
# Segment CA ends
# Segment DA begins
####################### Segment SC begins as image of DA #######################
set SC [atomselect $m1 "serial 12071 to 12075"]
$SC set segname SC
set SC_data [ backup $SC [ list chain x y z resid resname name ] ]
$SC set chain SC
$SC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$SC set resid [list 606 606 606 606 606]
$SC writepdb segtype_generic_SC.pdb
segment SC {
    first none
    last none
    pdb segtype_generic_SC.pdb
}
coordpdb segtype_generic_SC.pdb SC
######################### Restoring A.U. state for DA ##########################
restore $SC [ list chain x y z resid resname name ]  $SC_data
############################### Segment SC ends ################################
# Segment DA ends
# Segment EA begins
####################### Segment TC begins as image of EA #######################
set TC [atomselect $m1 "serial 12076 to 12080"]
$TC set segname TC
set TC_data [ backup $TC [ list chain x y z resid resname name ] ]
$TC set chain TC
$TC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$TC set resid [list 607 607 607 607 607]
$TC writepdb segtype_generic_TC.pdb
segment TC {
    first none
    last none
    pdb segtype_generic_TC.pdb
}
coordpdb segtype_generic_TC.pdb TC
######################### Restoring A.U. state for EA ##########################
restore $TC [ list chain x y z resid resname name ]  $TC_data
############################### Segment TC ends ################################
# Segment EA ends
# Segment IA begins
####################### Segment UC begins as image of IA #######################
set UC [atomselect $m1 "serial 12123 to 12127"]
$UC set segname UC
set UC_data [ backup $UC [ list chain x y z resid resname name ] ]
$UC set chain UC
$UC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$UC set resid [list 701 701 701 701 701]
$UC writepdb segtype_generic_UC.pdb
segment UC {
    first none
    last none
    pdb segtype_generic_UC.pdb
}
coordpdb segtype_generic_UC.pdb UC
######################### Restoring A.U. state for IA ##########################
restore $UC [ list chain x y z resid resname name ]  $UC_data
############################### Segment UC ends ################################
# Segment IA ends
# Segment JA begins
####################### Segment VC begins as image of JA #######################
set VC [atomselect $m1 "serial 12128 to 12132"]
$VC set segname VC
set VC_data [ backup $VC [ list chain x y z resid resname name ] ]
$VC set chain VC
$VC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$VC set resid [list 702 702 702 702 702]
$VC writepdb segtype_generic_VC.pdb
segment VC {
    first none
    last none
    pdb segtype_generic_VC.pdb
}
coordpdb segtype_generic_VC.pdb VC
######################### Restoring A.U. state for JA ##########################
restore $VC [ list chain x y z resid resname name ]  $VC_data
############################### Segment VC ends ################################
# Segment JA ends
# Segment MA begins
####################### Segment WC begins as image of MA #######################
set WC [atomselect $m1 "serial 12152 to 12168"]
$WC set segname WC
set WC_data [ backup $WC [ list chain x y z resid resname name ] ]
$WC set chain WC
$WC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$WC set resid [list 701 702 703 704 705 706 707 708 709 710 711 712 713 714 715 716 717]
$WC writepdb segtype_generic_WC.pdb
segment WC {
    first none
    last none
    pdb segtype_generic_WC.pdb
}
coordpdb segtype_generic_WC.pdb WC
######################### Restoring A.U. state for MA ##########################
restore $WC [ list chain x y z resid resname name ]  $WC_data
############################### Segment WC ends ################################
# Segment MA ends
# Segment NA begins
####################### Segment XC begins as image of NA #######################
set XC [atomselect $m1 "serial 12169 12170"]
$XC set segname XC
set XC_data [ backup $XC [ list chain x y z resid resname name ] ]
$XC set chain XC
$XC move { { -0.5 0.8660254038 0.0 -257.78  } { -0.8660254038 -0.5 0.0 -446.4880571751  } { 0.0 0.0 1.0 0.0  } { 0.0 0.0 0.0 1.0  }  }
$XC set resid [list 801 802]
$XC writepdb segtype_generic_XC.pdb
segment XC {
    first none
    last none
    pdb segtype_generic_XC.pdb
}
coordpdb segtype_generic_XC.pdb XC
######################### Restoring A.U. state for NA ##########################
restore $XC [ list chain x y z resid resname name ]  $XC_data
############################### Segment XC ends ################################
# Segment NA ends
############################# DISU patches follow ##############################
patch DISU E:24 E:44
patch DISU E:89 E:175
patch DISU E:96 E:166
patch DISU E:101 E:119
patch DISU E:188 E:217
patch DISU E:198 E:209
patch DISU E:266 E:300
patch DISU E:347 E:413
patch DISU E:354 E:386
patch DISU F:87 F:93
############################# LINK patches follow ##############################
patch NGLA E:58 EG01:1
patch NGLB E:103 EG18:1133
patch NGLA E:107 EG02:1
patch NGLA E:118 EG09:1
patch NGLA E:122 EG10:1
patch NGLA E:167 EG11:1
patch NGLA E:204 EG03:1
patch NGLA E:232 EG04:1
patch NGLA E:246 EG15:1276
patch NGLA E:265 EG12:1
patch NGLB E:271 EG13:1
patch NGLA E:301 EG14:1
patch NGLB E:308 EG16:1839
patch NGLA E:324 EG17:1355
patch NGLA E:332 EG05:1
patch NGLA E:355 EG06:1
patch NGLB E:361 EG07:1
patch NGLA E:416 EG08:1
patch NGLB F:100 FG01:1611
patch NGLA F:107 FG02:1618
patch NGLA F:126 FG03:1637
patch 14ab EG01:1 EG01:2
patch 14ab EG01:2 EG01:3
patch 16BT EG01:3 EG01:4
patch 13aa EG01:3 EG01:7
patch 13ba EG01:4 EG01:5
patch 16AT EG01:4 EG01:6
patch 14ab EG02:1 EG02:2
patch 14ab EG02:2 EG02:3
patch 13ba EG02:3 EG02:4
patch 14ab EG03:1 EG03:2
patch 14bb EG04:1 EG04:2
patch 14aa EG04:2 EG04:3
patch 13ba EG04:3 EG04:4
patch 16AT EG04:3 EG04:6
patch 12ba EG04:4 EG04:5
patch 14bb EG05:1 EG05:2
patch 14bb EG06:1 EG06:2
patch 14bb EG07:1 EG07:2
patch 14ab EG08:1 EG08:2
patch 14bb EG08:2 EG08:3
patch 14ab EG09:1 EG09:2
patch 14ab EG09:2 EG09:3
patch 13ba EG09:3 EG09:4
patch 16AT EG09:3 EG09:5
patch 14ab EG10:1 EG10:2
patch 14ab EG11:1 EG11:2
patch 14ab EG12:1 EG12:2
patch 14bb EG13:1 EG13:2
patch 14bb EG14:1 EG14:2
patch 14bb EG14:2 EG14:3
patch 13aa EG14:3 EG14:4
patch 16AT EG14:3 EG14:7
patch 12aa EG14:4 EG14:5
patch 12ab EG14:5 EG14:6
patch 13bb EG14:7 EG14:8
patch 16AT EG14:7 EG14:10
patch 12bb EG14:8 EG14:9
############################### Transform 2 ends ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap 00-01-00_psfgen-build.psf
writepdb 00-01-00_psfgen-build.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
