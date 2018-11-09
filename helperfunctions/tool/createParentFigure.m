function handle = createParentFigure(handle)

handle = figure('position', [100, 100, 1000, 500],'Name','pipeline');

a = findall(handle);
for i=1:length(a)
    if isa(a(i),'matlab.ui.container.Menu')
        set(a(i),'Visible','Off')
    end
end

b = findall(a,'Tag','Annotation.InsertLegend');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','Standard.NewFigure');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','Standard.FileOpen');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','DataManager.Linking');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','Annotation.InsertColorbar');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','Plottools.PlottoolsOn');
for i=1:length(b)
     set(b,'Visible','Off')
end

b=findall(a,'Tag','Plottools.PlottoolsOff');
for i=1:length(b)
     set(b,'Visible','Off')
end

end