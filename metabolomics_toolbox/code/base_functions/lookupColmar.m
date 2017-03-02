function names=lookupColmar(in,db)
%in is an array of doubles with shifts
%out is a cell array of strings
%ex: names=lookupColmar([4.02,3.58,3.49,3.23],'H')

%construct the url
url1='http://spinportal.magnet.fsu.edu/webquery/webquery.php?peaklist=[';
url2=']&numreturns=10&database=';
url3='&detail=no';

if exist('db','var') & upper(db)=='C'
    dburl='BMRB_C13';
elseif exist('db','var') & upper(db)=='H'
    dburl='BMRB';
else
    error('Need database, ''H'' or ''C''')
end

in=num2str(in,'%0.1f,');
in(end)=[];
in(in==' ')='';

url=[url1 in url2 dburl url3];

page=urlread(url);

%parse page
idx=strfind(page,'Input Peaklist');
page=strread(page(idx:idx+1000), '%s', 'delimiter', sprintf('\n'));
page(1)=[];

for i=1:length(page)
    a=strsplit(page{i});
    names{i}=a{1};
end

i=1;
while i<=length(names)
    if names{i}(1)=='<'
        names(i:end)=[];
        i=length(names);
    end
    i=i+1;
end
