function vals = Read_Tecmag_Header_Var(filename,Var_names)
% Modified from Read_Tecmag_Header.m by Dave Korenchan
% header=[];
    
if exist(filename,'file')
    s=dir(filename);
    if s.bytes>1056
%         Hdr_var=Read_Tecmag_hdr;
% 
%         Hdr_name=Hdr_var(:,1:28);
%         Hdr_offset=Hdr_var(:,29:34);
%         Hdr_type=Hdr_var(:,36:43);
%         Hdr_size=Hdr_var(:,45:47);
%         Hdr_desc=Hdr_var(:,48:85);
%         % function strrep pour enlever les espaces.
%     
        Number_var=numel(Var_names);
        vals=cell(Number_var,1);

        % Read entire content of file
        fcontent=readlines(filename);

        % Then search through file content for each instance of the
        % variable name. The pattern to find is: the variable name, then a
        % combo of whitespace/control characters, then some
        % non-whitespace/non-control combo of characters
        for ii=1:Number_var
            Var_name=Var_names{ii};
            pattern=[Var_name,'[\s\x00-\x1F]+','[^\s\x00-\x1F]+'];
            matches=regexp(fcontent,pattern,'match');
            matches=char([matches{:}].'); %convert from str to char
            matches(1,:)=[];  % we can throw out the first match; I think 
                            % this pertains to something basic about the
                            % sequence...

            % Now, for each identified match, ID if it's a mathematical
            % expression or a value, and save in either tempexprs or
            % tempvals, respectively
            tempexprs=cell(size(matches,1),1);
            tempvals=cell(size(matches,1),1);
            for jj=1:size(matches,1)
                entry=matches(jj,:);
                val=extractAfter(entry,Var_name);
                if contains(val,'=')
                    val=val(find(~isstrprop(val,'cntrl')&~isstrprop(val,'wspace')));
                    tempexprs{jj}=val;
                else
                    tv=Parse_Tecmag_Variable(val);
                    if ~isnan(tv)
                        tempvals{jj}=tv;
                    end
                end
            end

            % Parse through obtained values/expressions. Raise a warning if
            % both a value and expression were obtained, or if
            % values/expressions are not equivalent
            if sum(~cellfun('isempty',tempexprs))>0 && sum(~cellfun('isempty',tempvals))>0
                warning(['Variable ' Var_name ' has both a value and an expression specified! ' ...
                    'We will return the expression.'])
                exprflg=true;
                evalvals=tempexprs(~cellfun('isempty',tempexprs));
            elseif sum(~cellfun('isempty',tempexprs))>0
                exprflg=true;
                evalvals=tempexprs(~cellfun('isempty',tempexprs));
            elseif sum(~cellfun('isempty',tempvals))>0
                exprflg=false;
                evalvals=tempvals(~cellfun('isempty',tempvals));
            end

            if exprflg
                if numel(unique(evalvals))>1
                    warning(['Variable ' Var_name ' has multiple specified expressions! ' ...
                        'We will use the first one found.'])
                end           
            else
                if numel(unique(cell2mat(evalvals)))>1
                    warning(['Variable ' Var_name ' has multiple specified values! ' ...
                        'We will use the first one found.'])
                end           
            end
            vals{ii}=evalvals{1};
        end

        % Now go through the extracted values, ID which contain
        % expressions, and evaluate them. If they cannot be evaluated,
        % return a warning to the user
        exprIdx=find(cellfun(@ischar,vals));
        for ii=1:numel(exprIdx)
            expr=vals{exprIdx(ii)};
            depVarNames=extractBetween(expr,'[',']');
            if sum(strcmp(Var_names,depVarNames))<numel(depVarNames)
                depVarNamesStr=[depVarNames{:}];
                warning(['Expression found for variable ' Var_names{exprIdx(ii)} ...
                    ' depends on variable(s): ' depVarNamesStr '! Cannot evaluate; '...
                    'returning expression string instead...'])
            else
                for jj=1:numel(depVarNames)
                    depVarNm=depVarNames{jj};
                    depVarVal=vals{strcmp(Var_names,depVarNm)};
                    expr=replace(expr,['[' depVarNm ']'],num2str(depVarVal));
                end
                expr=extractAfter(expr,'='); %eval() doesn't like = signs
                vals{exprIdx(ii)}=eval(expr);
            end
        end
        
%         File_pt=fopen(filename,'r','l');
%         remCarryOver=false;
%         rem=[];
%         while 1
%             if remCarryOver
%                 line=rem;
%                 remCarryOver=false;
%             else
%                 line=fgetl(File_pt);
%             end
%             for i=1:Number_var
%                 Var_name=Var_names{i};
%                 if strfind(line,Var_name)>0
%                     %Make sure character preceding Var_name isn't an =
%                     %sign! That means it's being called for something else
%                     if ~strcmp(line(strfind(line,Var_name)),'=')
%                         %Next non-whitespace, non-control character after 
%                         %Var_name should be either the value, or an 
%                         %expression (starting with an = sign)
%                         rem=line(strfind(line,Var_name)+length(Var_name):end);
%                         nextrem=rem(find(~isstrprop(rem,'cntrl')&~isstrprop(rem,'wspace')):end);
%                         val=nextrem(1:find(isstrprop(nextrem,'cntrl')|isstrprop(nextrem,'wspace'))-1);
%                         if isstrprop(val(1),'digit') %obtain + parse value
%                             vals{i}=Parse_Tecmag_Variable(val);
%                         elseif strcmp(val(1),'=') %obtain expression; evaluate after all variables read in
%                             expr{i}=val;
%                         end
%                     end
%                 end
%             end
%             if ~ischar(line) || isempty(Var_names)
%                 break
%             end
                
                
%                 Var_offset=eval(Hdr_offset(i,:));
%                 Var_type=strrep(Hdr_type(i,:),' ','');
%                 Var_size=eval(Hdr_size(i,:));
%                 Var_desc=Hdr_desc(i,:);
% 
%                 fseek(File_pt,Var_offset,'bof');
%                 Var_data=fread(File_pt,Var_size,Var_type)';
    %             if (strcmp(Var_type,'char'))&((Var_size-1)>0)
    %                 Var_data=char(Var_data);
    %             end
%                 eval(strcat('header.',Var_name,'=Var_data;'));
    end
%     fclose(File_pt);
end
end