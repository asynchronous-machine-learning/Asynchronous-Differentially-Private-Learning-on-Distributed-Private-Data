clc
clear

%#ok<*SAGROW>
%#ok<*NASGU>
%#ok<*VUNUS>
%#ok<*NOPTS>
%#ok<*FNDSB>

fileID     =fopen('hospital.csv');
fileID_outX=fopen('X.csv','w');
fileID_outY=fopen('Y.csv','w');
%%

fgetl(fileID);
line=1;
i=1;
while(~feof(fileID))
    i=i+1;
    line=fgetl(fileID);
    place_quation=find(line=='"');
    while(length(place_quation)>1)
        line=strrep(line,line(place_quation(1):place_quation(2)),' ');
        place_quation=find(line=='"');
    end
    parsed_line=split(line,',');

    table(1,1)=str2double(parsed_line{4});    
    table(1,2)=str2double(parsed_line{6}(1:2));
    
    switch parsed_line{8}
        case 'M'
            table(1,3)=1;
        case 'F'
            table(1,3)=2;
        otherwise
            table(1,4)=3;
    end
    
    switch parsed_line{9}
        case 'Black/African American'
            table(1,4)=1;
        case 'Multi'
            table(1,4)=2;
        case 'Other Race'
            table(1,4)=3;
        case 'Unknown'   
            table(1,4)=4;
        case 'White'
            table(1,4)=5;
        otherwise
            table(1,4)=6;
    end    

    switch parsed_line{10}
        case 'Spanish/Hispanic'
            table(1,5)=1;
        case 'Not Span/Hispanic'
            table(1,5)=2;
        case 'Multi'
            table(1,5)=3;
        case 'Unknown'   
            table(1,5)=4;
        otherwise
            table(1,5)=5;
    end  
    
    table(1,6)=str2double(parsed_line{11});
    
    switch parsed_line{12}
        case 'Emergency'
            table(1,7)=1;
        case 'Urgent'
            table(1,7)=2;
        case 'Elective'
            table(1,7)=3;
        case 'Newborn'   
            table(1,7)=4;
        case 'Not Available'
            table(1,7)=5;
        case 'Trauma'
            table(1,7)=6;
        otherwise
            table(1,7)=7;
    end      
    
    table(1,8)=str2double(parsed_line{15});
    
    table(1,9)=str2double(parsed_line{17});
    
    table(1,10)=str2double(parsed_line{19});
    
    table(1,11)=str2double(parsed_line{21});
    
    table(1,12)=str2double(parsed_line{23});
    
    Y(1,1)=str2double(parsed_line{36});
    
    fprintf(fileID_outX,'%s\n',regexprep(num2str(table),'\s+',','));
    fprintf(fileID_outY,'%s\n',regexprep(num2str(Y)    ,'\s+',','));
    if i/1000==floor(i/1000)
        i
    end
    if sum(abs(imag(table)))
        break
    end
end