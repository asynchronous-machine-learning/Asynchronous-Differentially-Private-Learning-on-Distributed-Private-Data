clc
clear

%#ok<*SAGROW>
%#ok<*NASGU>
%#ok<*VUNUS>
%#ok<*NOPTS>
%#ok<*FNDSB>

fileID     =fopen('loan.csv');
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
    while(length(parsed_line)>74)
        for j=21:length(parsed_line)-1
            parsed_line{j}=parsed_line{j+1};
        end
        parsed_line(end)=[];
    end
    for j=3:length(parsed_line)
        if j<6
            table(1,j-2)=str2double(parsed_line{j});
            if isnan(table(1,j-2))
                table(1,j-2)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end
        if j==6
            parsed_line{j}=erase(parsed_line{j},'months');
            table(1,j-3)=str2double(parsed_line{j}(find(~isspace(parsed_line{j}))));
        end
        if j==7
            Y(1,1)=str2double(parsed_line{j});
        end
        if (j>7) && (j<9)
            table(1,j-3)=str2double(parsed_line{j});
            if isnan(table(1,j-3))
                table(1,j-3)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end
        if j==9
            table(1,j-3)=double(parsed_line{j})-double('A');
        end
        if j==10
            table(1,j-3)=str2double(erase(parsed_line{j},parsed_line{j-1}));
        end
        if j==12
            parsed_line{j}=erase(parsed_line{j},'+ years');
            parsed_line{j}=erase(parsed_line{j},'years');
            parsed_line{j}=erase(parsed_line{j},'year');
            parsed_line{j}=parsed_line{j}(find(~isspace(parsed_line{j})));
            if ~strcmp(parsed_line{j},'n/a')
                if strcmp(parsed_line{j}(1),'<')
                    table(1,j-4)=0;
                else
                    table(1,j-4)=str2double(parsed_line{j});
                end
            else
                table(1,j-4)=-1;
            end
        end 
        if j==13
            switch parsed_line{j}
                case 'RENT'
                    table(1,j-4)=1;
                case 'OWN'
                    table(1,j-4)=2;
                case 'MORTGAGE'
                    table(1,j-4)=3;
                otherwise
                    table(1,j-4)=4;
            end
        end
        if j==14
            table(1,j-4)=str2double(parsed_line{j});
        end       
        if j==15
            switch parsed_line{j}
                case 'Verified'
                    table(1,j-4)=1;
                case 'Source Verified'
                    table(1,j-4)=2;
                case 'Not Verified'
                    table(1,j-4)=3;
                otherwise
                    table(1,j-4)=4;
            end
        end
        if j==16
            if ~isempty(parsed_line{j})
                table(1,j-4)=datenum(parsed_line{j},'mmm-yyyy')-datenum('Feb-1970','mmm-yyyy');
            else
                table(1,j-4)=0;
            end            
        end        
        if j==17
            switch parsed_line{j}
                case 'Fully Paid'
                    table(1,j-4)=1;
                case 'Charged Off'
                    table(1,j-4)=2;
                case 'Current'
                    table(1,j-4)=3;
                otherwise
                    table(1,j-4)=4;
            end            
        end
        if j==18
            table(1,j-4)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j})))))-double('n');
        end
        if (j>21) && (j<23)
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end        
        if j==23
            parsed_line{j}=erase(parsed_line{j},' ');
            parsed_line{j}=erase(parsed_line{j},'xx');
            table(1,j-6)=str2double(parsed_line{j});
        end   
        if (j>23) && (j<27)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end   
        if j==27
            if ~isempty(parsed_line{j})
                table(1,j-6)=datenum(parsed_line{j},'mmm-yyyy')-datenum('Feb-1970','mmm-yyyy');
            else
                table(1,j-6)=0;
            end            
        end
        if (j>27) && (j<36)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end          
        if j==36
            table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j})))))-double('f');
        end
        if (j>36) && (j<46)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end          
        if j==46
            if ~isempty(parsed_line{j})
                table(1,j-6)=datenum(parsed_line{j},'mmm-yyyy')-datenum('Feb-1970','mmm-yyyy');
            else
                table(1,j-6)=0;
            end            
        end
        if (j==47)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end          
        if j==48
            if ~isempty(parsed_line{j})
                table(1,j-6)=datenum(parsed_line{j},'mmm-yyyy')-datenum('Feb-1970','mmm-yyyy');
            else
                table(1,j-6)=0;
            end
        end
        if j==49
            if ~isempty(parsed_line{j})
                table(1,j-6)=datenum(parsed_line{j},'mmm-yyyy')-datenum('Feb-1970','mmm-yyyy');
            else
                table(1,j-6)=0;
            end
        end
        if (j>49) && (j<53)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end          
        if j==53
            switch parsed_line{j}
                case 'INDIVIDUAL'
                    table(1,j-6)=1;
                otherwise
                    table(1,j-6)=2;
            end     
        end     
        if (j>53)
            parsed_line{j}=erase(parsed_line{j},' ');
            table(1,j-6)=str2double(parsed_line{j});
            if isnan(table(1,j-6))
                table(1,j-6)=sum(double(parsed_line{j}(find(~isspace(parsed_line{j}))))); 
            end
        end          
    end
    fprintf(fileID_outX,'%s\n',regexprep(num2str(table),'\s+',','));
    fprintf(fileID_outY,'%s\n',regexprep(num2str(Y),'\s+',','));
    if i/1000==floor(i/1000)
        i
    end
    if sum(abs(imag(table)))
        break
    end
end