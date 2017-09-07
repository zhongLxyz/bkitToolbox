function Response = recodeResponse(response)

% this function recodes the response column read from the data files and
% output a colunm with 1s and 0s.

% for responses coded in '0i' and 'i0', '0i' is recoded as 1, which means
% the response corresponding to the left arrow key is 1.
% (the order of the options on the screen: i0, 0i)
% (CW, CCW) (above, below) (towards, away)

if ismember('CCW', response) || ismember('CCW<br/>', response)
    Response = strncmpi('CC',response,2);
elseif ismember('''CCW''', response)
    Response = strncmpi('''CCW''',response,6);
elseif ismember('CounterClockwise', response) || ismember('Counterclockwise<br/>', response)
    Response = strncmpi('Counter',response,6);
elseif ismember('Left', response) || ismember('Left<br/>', response)
    Response = strncmpi('Le',response,2); % left = 1, right = 0
elseif ismember('0i',response)|| ismember('i0',response) || ismember('i0<br/>',response) || ismember('0i<br/>',response)
    Response = strncmpi('0i',response,2);
elseif ismember('''Away''',response) || ismember('''Toward''',response)
    Response = strncmpi('''Away''',response,4);
elseif ismember('''Above''',response)
    if ismember('''Below''',response)
        Response = strncmpi('''Below''',response,4);
    else
        Response = not(strncmpi('''Above''',response,4));
    end
elseif ismember('Above''',response)
    if ismember('Below''',response)
        Response = strncmpi('Below''',response,4);
    else
        Response = not(strncmpi('Above''',response,4));
    end
end

