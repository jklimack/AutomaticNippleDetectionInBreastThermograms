

function overlay = nipple_klimack(imgo, display)
    
    shape = size(imgo);
    if length(shape) == 3
        img = rgb2gray(imgo);
    end
    
    % STEP 1: HUMAN BODY SEGMENTATION
    % thresholding to create the body mask
    btheta = 50;
    bmask = img < btheta;
    % morphologically close mask using disc-shaped element with r=3
    disc = strel('disk', 3);
    bmask = imclose(bmask, disc);
    % morphological dilation with disc-shaped element with r=10
    disc = strel('disk', 10);
    bmask = imdilate(bmask, disc);
    % Invert thresholded image
    bmask = not(bmask);
    bseg = zeros(shape(1), shape(2), 'uint8');
    bseg(bmask) = img(bmask);
    bseg = double(bseg);
    bseg = bseg / 255.0;

    % STEP 2: NIPPLE CANDIDATES
    % perform adaptive thresholding
    % neighbourhood size is 15 pixels, and the statistic is median. 
    num_neighbourhood = 15;
    img_c = medfilt2(bseg, [num_neighbourhood, num_neighbourhood]);
    % subtract original image from convolved one and threshold the image
    x = img_c - bseg;
    C = 0.03;
    x = x > C;

    % STEP 3: SELECTION ALGORITHM
    % find regions for each candidate
    R = compute_regions(x); 
    % remove candidates in the upper and lower regions
    Rup = 0.35 * shape(1);
    Rlw = shape(1) - 0.3 * shape(1);
    Rext = zeros(shape(1),shape(2), 'logical');
    Rext(1:Rup,:) = 1;
    Rext(Rlw:shape(1),:) = 1;
    N_up_lw = zeros(shape(1), shape(2), 'int32');
    N_up_lw(Rext) = R(Rext);
    N_up_lw = unique(N_up_lw);
    for r = 1:size(N_up_lw)
       idx = find(R==N_up_lw(r));
       R(idx) = 0;
    end
    % remove candidates outside the body mask
    N_hb(bmask) = R(bmask);
    N_hb = unique(N_hb);
    for r = 1:size(N_hb)
       idx = find(R==N_hb(r));
       R(idx) = 0;
    end
    % split candidates into left and right regions
    Lcnt = shape(2) / 2;
    N_left = zeros(shape(1), shape(2), 'int32');
    N_right = zeros(shape(1), shape(2), 'int32');
    N_left(:,1:Lcnt) = R(:,1:Lcnt);
    N_right(:,Lcnt:shape(2)) = R(:,Lcnt:shape(2));
    % compute the left and right nipple targets
    RT = compute_nipple(N_right);
    LT = compute_nipple(N_left);
    
    
    % OVERLAY POINTS ON IMAGE
    if display
        image(imgo);
        hold on
        axis off;
        %plot([1 1]*Lcnt, ylim, '-y') % center line
        %plot(xlim, [1 1]*Rup, '-b') % upper region line
        %plot(xlim, [1 1]*Rlw, '-g') % lower region line
        scatter([RT(2) LT(2)], [RT(1) LT(1)], 'r*');
        hold off
    end
    
    overlay = {imgo; [RT(2) LT(2)]; [RT(1) LT(1)]};
end

function T = compute_nipple(Nm)
    shape = size(Nm);
    Nm = drop_small_regions(Nm);
    N = unique(Nm);
    N = N(2:size(N));
    if size(N) > 1
        % compute roundness of each region in N
        A = zeros(size(N), 'uint32');
        for i = 1:size(N)
            % extract region information
            region = zeros(shape(1), shape(2), 'int32');
            [ri, ci] = find(Nm==N(i));
            region([ri, ci]) = Nm([ri, ci]);
            region = int32(region > 0);
            Atemp = size(ri);
            A(i) = Atemp(1);
            % compute perimeter
            P = zeros(1);
            b = boundary(ri, ci);
            r = ri(b(end));
            c = ci(b(end));
            for j=1:size(b)
                P(1) = P(1) + sqrt((r-ri(b(j))).^2 + (c-ci(b(j))).^2);
                r = ri(b(j));
                c = ci(b(j));
            end
            % compute roundness
            R(i) = 4*pi*A(i) / (P(1)^2);
        end
        m = find(R==max(R))
        sm = size(m)
        if sm(2) >1
            maxA = find(A==max(A(m)));
            T = N(maxA(1));
        else
            T = N(m);
        end
    else
        s = size(N);
        if s(2) == 0
            T = -1;
        else
            T = N(1);
        end
    end
    if T == -1
        rm=0;cm=0;
    else
        [ri, ci] = find(Nm==T);
        rm = uint32(mean(ri));
        cm = uint32(mean(ci));
    end
    T = [rm, cm];
end

function regions = compute_regions(x)
    shape = size(x);
    num_pos = sum(x, 'all');
    [rows,cols] = find(x);
    x = int32(x);
    
    R = 1;
    for i=1:size(rows)
       if x(rows(i),cols(i))==1 | x(rows(i),cols(i))==-1
           R = R+1;
           Sr = 0;Sc=0;
           clear Sr Sc
           Sr(1) = rows(i);
           Sc(1) = cols(i);
           p1 = 1;
           p2 = 2;
           while not(p1==p2)
               r = Sr(p1);
               c = Sc(p1);
               p1 = p1+1;
               x(r,c) = R;
               
               if not(r==shape(1))
                   if x(r+1,c)==1 %right
                      x(r+1,c)=-1;
                      Sr(p2) = r+1;
                      Sc(p2) = c;
                      p2 = p2+1;
                   end
               end
               if not(r==1)
                   if x(r-1,c)==1 %left
                      x(r-1,c)=-1;
                      Sr(p2) = r-1;
                      Sc(p2) = c;
                      p2 = p2+1;
                   end
               end
               if not(c==1)
                   if x(r,c-1)==1 %up
                      x(r,c-1)=-1;
                      Sr(p2) = r;
                      Sc(p2) = c-1;
                      p2 = p2+1;
                   end
               end
               if not(c==shape(2))
                   if x(r,c+1)==1 %down
                      x(r,c+1)=-1;
                      Sr(p2) = r;
                      Sc(p2) = c+1;
                      p2 = p2+1;
                   end
               end
           end
       end
    end
    regions = x;
end 

function R = drop_small_regions(x)
    N = unique(x);
    N = N(2:size(N));
    
    % sort by area size in decreasing order
    sn = size(N);
    A = zeros(sn(2), 2,'uint32');
    for r = 1:size(N)
       idx = find(x==N(r));
       si = size(idx);
       A(r,1) = r;
       A(r,2) = si(1);
    end
    A = sortrows(A, 2, 'descend');
    
    % remove candidates smaller than Np=20
    Np = 20;
    for i = 1:size(N)
       if A(i,2) < Np & i>1
           idx = find(x==N(A(i,1)));
           x(idx) = 0;
       end
    end
    R = x;
end

