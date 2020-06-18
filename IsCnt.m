function IsCnt

disp("[Press: 'd' to continue, 'space' to break]");
stp = true;
while stp == true
    w = waitforbuttonpress; 
    switch w 
        case 1 % (keyboard press) 
          key = get(gcf,'currentcharacter'); 
              switch key
                  case 100 % 100 is lowercase d
                      disp('Proceeding ...')
                      stp = false;
                  case 32 % 32 is the space key
                      disp('Breaking the Analysis ...')
                      break % break out of the while loop
                  otherwise
                      disp('Press a Correct Key!')
             end
     end
end
 
end

