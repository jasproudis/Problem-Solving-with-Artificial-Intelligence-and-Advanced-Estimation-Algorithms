%eisagwgh iris dataset:
fid = fopen('iris.data');  #prosoxi edw vazoume to absolute path tou iris.data ston diko mas upologisti
iris_data = textscan(fid, '%f %f %f %f %s', 200, 'Delimiter', ',');
fclose(fid);
irisMatrix = cell2mat(iris_data(:,1:4));

%afairoume thn teleutaia grammh se periptwsh pou einai NaN
irisMatrix(end,:) = [];
disp(irisMatrix);

%statistikoi upologismoi (mesh timh, tupikh apoklish)

mu_vector = mean(irisMatrix); %epistrefei ena dianusma me tis 4 meses times twn sthlwn
printf("\n");
%disp(mu_vector);

%apo8hkeuoume kathe mesh timh se ksexwristh metavlhth
mu_sl = mu_vector(1);
mu_sw = mu_vector(2);
mu_pl = mu_vector(3);
mu_pw = mu_vector(4);
printf("Sepal Lenght mean = %d\n", mu_sl);
printf("Sepal Width mean = %d\n", mu_sw);
printf("Petal Lenght mean = %d\n", mu_pl);
printf("Petal Width mean = %d\n", mu_pw);

std_vector = std(irisMatrix); %epistrefei ena dianusma me 4 tupikes apokliseis (twn sthlwn)
disp(std_vector);
printf("\n");
%apo8hkeuoume kathe tupikh apoklish se ksexwristh metavlhth
std_sl = std_vector(1);
std_sw = std_vector(2);
std_pl = std_vector(3);
std_pw = std_vector(4);
printf("Sepal Lenght std = %d\n", std_sl);
printf("Sepal Width std = %d\n", std_sw);
printf("Petal Lenght std = %d\n", std_pl);
printf("Petal Width std = %d\n", std_pw);


%h subplot mas epitrepei na exoume polles diaforetikes sunarthseis ektupwmenes se ena figure
%typwnoume ta istogrammata twn 4 sthlwn tou dataset
subplot(4,1,1)
hist(irisMatrix(:,1), 15);
subplot(4,1,2)
hist(irisMatrix(:,2), 15);
subplot(4,1,3);
hist(irisMatrix(:,3), 15);
subplot(4,1,4);
hist(irisMatrix(:,4), 15);

figure(2);

%kanoyme fit tis 4 sthles twn metavlhtwn se mia kanonikh katanomh Gauss

subplot(4,1,1)
f1 = exp(-(irisMatrix(:,1)-mu_sl).^2./(2*std_sl^2))./(std_sl*sqrt(2*pi));
plot(irisMatrix(:,1), f1, 'o');

subplot(4,1,2)
f2 = exp(-(irisMatrix(:,2)-mu_sw).^2./(2*std_sw^2))./(std_sw*sqrt(2*pi));
plot(irisMatrix(:,2), f2, 'o');

subplot(4,1,3)
f3 = exp(-(irisMatrix(:,3)-mu_pl).^2./(2*std_pl^2))./(std_pl*sqrt(2*pi));
plot(irisMatrix(:,3), f3, 'o');

subplot(4,1,4)
f4 = exp(-(irisMatrix(:,4)-mu_pw).^2./(2*std_pw^2))./(std_pw*sqrt(2*pi));
plot(irisMatrix(:,4), f4, 'o');


%dhmiourgoume ton pinaka kai ton suntelesth susxetishs
figure(3);
plotmatrix(irisMatrix);
rho = corr(irisMatrix);
disp(rho);

