function writeDatas(t, Wc, Vc, Wb, Vb, r, dNt)
    data(1)=libpointer('doublePtrPtr', t);
    data(2)=libpointer('doublePtrPtr', Wc);
    data(3)=libpointer('doublePtrPtr', Vc);
    data(4)=libpointer('doublePtrPtr', Wb);
    data(5)=libpointer('doublePtrPtr', Vb);

    files = {};
    files{1} = ['����\t', num2str(r), '.csv'];
    files{2} = ['����\Wc', num2str(r), '.csv'];
    files{3} = ['����\Vc', num2str(r), '.csv'];
    files{4} = ['����\Wb', num2str(r), '.csv'];
    files{5} = ['����\Vb', num2str(r), '.csv'];
    
    for j = 1: length(data)
        d = get(data(j), 'Value')';
        dlmwrite(files{j}, d(1: dNt: end, :), 'delimiter', ',', 'precision', 8, 'newline', 'pc');
    end
end