#CREATE STARTING BOARD
pieces <- c("*","*","*","*","*","*","*","*","*") 
board <- matrix(pieces,nrow = 3, ncol = 3, byrow = TRUE)
print(board)

#--------------------------------------------
#code to choose x or o
repeat {
  
  msg <-"\nX or O? (select 1 for x, 2 for o)"
  
  if(interactive()) {
    choice <- readline(msg)
  } else {
    cat(msg)
    choice <- readLines("stdin", 1)
  }
  
  if (choice %in% c(1,2) ) {
    if (choice == 1) {
      piece <- "x"
      piececomp <- "o" #piececomp used for computer piece
      cat("you have choosen x\n")
    } else if (choice == 2) {
      piece <- "o"
      piececomp <- "x"
      cat("you have choosen o\n")
    }
    break 
  } else if (!(choice %in% c(1,2)) | choice == "noerror") {
    cat("\nwrong choice, choose only 1 for x or 2 for o\n")
  } 
}

#-------------------------------------------------------------------------------
#code to play a turn

play <- function (posrow="noerror",poscol="noerror",...) {
  
  repeat {
    #allows code to loop through play if it is not a stalemate
    if (sum(board == "*") > 1) {
      #-------------------------------------------------------------------------
      #interactive code to choose row
      msgrow <-"\nChoose the row you would like your piece to be on. Only choose options 1 to 3\n"
      
      if(interactive()) {
        posrow <- as.integer(readline(msgrow))
      } else {
        cat(msgrow)
        posrow <- as.integer(readLines("stdin", 1))
      }
      
      cat("the row you chose is", posrow, "\n")
      
      #interactive code to choose column
      msgcol <-"\nChoose the column you would like your piece to be on. Only choose options 1 to 3\n"
      
      if(interactive()) {
        poscol <- as.integer(readline(msgcol))
      } else {
        cat(msgcol)
        poscol <- as.integer(readLines("stdin", 1))
      }
      
      cat("the column you chose is", poscol, "\n")
      
      suppressWarnings(as.numeric(posrow))
      suppressWarnings(as.numeric(poscol))
      #-------------------------------------------------------------------------
      if (posrow == "noerror" | poscol == "noerror" | !(posrow %in% c(1,2,3)) | !(poscol %in% c(1,2,3))) {
        print("wrong input. 1) only input values 1,2 or 3 for either row or column to match a coordinate on the board")
        next 
      } else {
        for (row in 1:nrow(board)) {
          for (col in 1:ncol(board)) {
            if (board[posrow,poscol] != "*") {
              print("this position is already taken, make another choice")
            } else if (board[posrow,poscol] == "*") {
              #player turn
              newboard <- board #make new board with changes that player makes
              newboard[posrow,poscol] <- piece #add player piece to position chosen
              
              #comp turn
              rows <- c(1,2,3)
              cols <- c(1,2,3)
              comprow <- sample(rows,1) #randomly assign computer position
              compcol <- sample(cols,1)
              if (newboard[comprow,compcol] != "*") {
                repeat { #ensures that the computer can select an open spot 
                  comprow <- sample(rows,1) #randomly assign computer position
                  compcol <- sample(cols,1)
                  if (newboard[comprow,compcol] == "*") {
                    newboard[comprow,compcol] <- piececomp
                    board <- newboard #rewrite board with newboard outside of function
                    print(board)
                    break
                  }
                }
              } else if (newboard[comprow,compcol] == "*") {
                newboard[comprow,compcol] <- piececomp
                board <- newboard #rewrite board with newboard outside of function
                print(board)
              }
            }
            break #breaks the for col loop 
          }
          #---------------------------------------------------------------------
          #code to decide winner/stalemate
          
          #checks if winner in rows
          win <- "no" #keeps track of wins to ensure that stalemate message shown only if there isn't a winner
          loss <- "no" #keeps track of losses so that code can be stopped if there is a loss before a stalemate
          for (wrow in 1:nrow(board)) {
            if (sum(board[wrow,]==piece) == 3) { 
              print ("You Won! Congrats!")
              win <- "yes"
              break 
            } else if (sum(board[wrow,]==piececomp) == 3) {
              win <- "no"
              loss <- "yes"
              print ("You Lost, Better Luck Next Time.")
            } 
          }
          
          #checks if winner in cols
          for (wcol in 1:ncol(board)) {
            if (sum(board[,wcol]==piece) == 3) {
              win <- "yes"
              print ("You Won! Congrats!")
              break 
            } else if (sum(board[,wcol]==piececomp) == 3) {
              win <- "no"
              loss <- "yes"
              print ("You Lost, Better Luck Next Time.")
            } 
          }
          
          #check if diagonal (L-R) win or loss
          if (board[1,1] == piece & board[2,2] == piece & board[3,3] == piece) {
            win <- "yes"
            print ("You Won! Congrats!")
          } else if (board[1,1] == piececomp & board[2,2] == piececomp & board[3,3] == piececomp) {
            win <- "no"
            loss <- "yes"
            print ("You Lost, Better Luck Next Time.")
          }
          
          #check if diagonal (R-L) win or loss         
          if (board[1,3] == piece & board[2,2] == piece & board[3,1] == piece) {
            win <- "yes"
            print ("You Won! Congrats!")
          } else if (board[1,3] == piececomp & board[2,2] == piececomp & board[3,1] == piececomp) {
            win <- "no"
            loss <- "yes"
            print ("You Lost, Better Luck Next Time.")
          }
          
          #check if stalemate
          if (win != "yes" & loss != "yes" & sum(board =="*") == 1) { #9 positions used with no winner yet (stalemate)
            print("STALEMATE -> no more available positions, the game ends in a tie")
            suppressWarnings(as.numeric(posrow))
            suppressWarnings(as.numeric(poscol))
          }
          #---------------------------------------------------------------------
          break #breaks the for row loop
        }
      }
    } else {
      break
    }
    #breaks code if there is a win regardless of whether or not its a stalemate
    if (win == "yes") {
      break
      #breaks code if there is a loss without a stalemate
      #win and sum board aren't in the same line bc play needs to be called for the correct win to be used 
    } else if (loss == "yes" & sum(board =="*") > 1) {
      break
    }
  }
}

play (posrow,poscol)