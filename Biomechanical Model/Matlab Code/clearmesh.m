function [node, elements] = clearmesh( node, elements )
t1 = tic;
    a=(1:size(node,1));
    empty_list = [];
    for i=1:size(node,1)
        if isempty(find(elements==a(i),1))
            disp(strcat('is empty :', num2str(i)))
            empty_list = [empty_list, i];
        end
    end
    
    empty_list = sort(empty_list,'descend');
    node = node(unique(elements),:);
    for i=1:numel(empty_list)
        disp( num2str( empty_list(i) ) );
       for j=1:numel(elements)
           if (elements(j))>(empty_list(i))
               elements(j) = elements(j)-1;
           end
       end
    end
disp(strcat('Time: :', num2str(toc( t1 ))));