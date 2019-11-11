<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
    <style>
        div.padded {
            padding-top: 0px;
            padding-right: 100px;
            padding-bottom: 0.25in;
            padding-left: 100px;
        }
    </style>
    <title>Final Project Report | CS 184</title>
    <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
    <link rel="stylesheet" type="text/css" href="style.css" media="screen"/>
    <script type="text/javascript"
            src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>
    <script type="text/x-mathjax-config">
  MathJax.Hub.Config({
    extensions: ["tex2jax.js"],
    jax: ["input/TeX", "output/HTML-CSS"],
    tex2jax: {
      inlineMath: [ ['$','$'], ["\\(","\\)"] ],
      displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
      processEscapes: true
    },
    "HTML-CSS": { availableFonts: ["TeX"] }
  });
</script>
    <script type="text/javascript" src="path-to-MathJax/MathJax.js">
    </script>
</head>
<body>
<br/>

<!-- Title, Summary and Team Members -->
<h1 align="middle">Final Project Report - Light Field Camera</h1>
<br>
<div class="padded" style="font-size:120%">
    <h3 align="middle">Members</h3>
    <p align="middle">Name: JIYUAN LU &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp SID: 3033483733</p>
    <p align="middle">Name: RENJIE SHAO &nbsp &nbsp &nbsp &nbsp SID: 3033530192</p>
    <p align="middle">Name: NAN WEI &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp SID: 3033530205</p>

    <h3 align="middle">Abstract</h3>
    <p>
        During project 3, we implemented a renderer that can trace radiance to generate the 2D image. During this process, light from the environment or objects is received by the lens and drops into pixels on sensor. The traditional way is to store the average spectrum received by each pixel. However, directional information of light is lost during this process. In order to reconstruct the light field, which can be represented by the function from the direction of light to the spectrum, we change each pixel to a grid that stores the direction information.
    </p>

    <p>
        With this extra information, our renderer can do refocusing, depth of field adjustment and camera position adjustment after rendering. All of these process can be finished in a few seconds, much faster than rendering. This technic definitely makes our images more editable.
    </p>

    <h3 align="middle">Technical Approach</h3>

    <h4>1. Store light field information</h4>

    <p>
        To store light field information, we change each pixel to a grid that stores the direction information.
    </p>

    <div align="center">
            <img src="images/light_field.png" width="600px"/>
            <figcaption align="middle">Real Light Field Camera</figcaption>
    </div>

    <p>
        In our model, microlens array is removed. Instead, we sample fixed positions on lens and store the radiance between each lens position and pixel pair.
    </p>

    <p>
        As is shown in the picture below, in order to store the light field, we need to know 4 position variables and the corresponding spectrum. <code>(u, v)</code> is the lens sample position and <code>(s, t)</code> is the pixel grid position. So we can reconstruct the light fild easily from this 4D representation.
    </p>

    <div align="center">
            <img src="images/Light_field_1.png" width="600px"/>
            <figcaption align="middle">Light Field Function</figcaption>
    </div>

    <h4>2. Refocusing</h4>
    <h4>2.1 Digital refocusing</h4>
    <p>
        <div align="center">
            <img src="images/refocus.png" width="400px"/>
            <figcaption align="middle">Fig 1</figcaption>
        </div>
        With the information of 4D light field, we can refocus images without doing ray tracing again. Given light field L(x,y,u,v), the irradiance on position (x,y) can be
        expressed as an integral:
        $$E(x,y) = \frac{1}{F^2}\int \int L(x,y,u,v)dudv$$
        $F$ is the distance between sensor plane and lens plane.
        The method to refocus is to change the separation between sensor plane and lens plane. As is depicted in Fig 1, change the separation from $F$ to $F'$ and the new light
        field between lens plane and new film plane can be expressed as:
        $$L'(x',y',u,v)=L(u+\frac{x'-u}{\alpha},v+\frac{y'-v}{\alpha}, u, v)$$
        $\alpha$ is $\frac{F'}{F}$. Therefore, for a position $(x',y')$ on new film plane, the irradiance can be expressed as:
        $$E'(x',y') = \frac{1}{\alpha^2F^2}\int\int L(u(1-\frac{1}{\alpha})+\frac{x'}{\alpha},
            v(1-\frac{1}{\alpha})+\frac{y'}{\alpha},u,v)$$
        when raytracing, we take $7*7$ samples on the lens for each pixel. To compute the refocus image, use the equation
        derived above and bilinear interpolation to compute irradiance
        at position $(x',y')$.
    </p>
    <h4>2.2 Anti-aliasing</h4>
    <p>
        Computed images will show step-edge aliasing in out of focus area which is caused
        by low sampling rate on lens. When the sampling rate is low,
        the shift in neighbouring subaperture images will larger than one pixel, causing step-edge artifacts.
        To anti-alias, we super sample the lens plane
        so that the images show less aliasing.
    </p>
    <h4>2.3 Depth estimation</h4>
    <p>
        For a better user experience, we try to allow users click on the
        image and obtain an image focusing on the selected point. To implement
        this, we have to estimate the depth of each position so that we know
        which focal length to apply when refocusing. The main idea is that if
        one position is out of focus, rays coming from different direction to
        this position will have a large variance. $\sigma_d(x,y)$ denotes variance of rays
        on position $(x,y)$ when focal plane is of depth $d$. The depth at each position
        (x,y) can be estimated as :
        $$depth(x,y)=argmax_d \sigma_d(x,y)$$
        We pre-processing light field to compute refocusing images for different focal plane and
        estimate depth for each position. To improve correctness of estimation, use a bilateral filter on $depth(x,y)$ to
        reduce noise.
    </p>
    <h4>3. Depth of field adjustment</h4>

    <p>
        To adjust depth of field, we use a simple method called <i>digitally stopping down the lens</i>. To get deeper depth of field, we decrease the number of samples we use to generate the final image. In effect, we get smaller aperture. In our experiment, we take 7*7 samples on lens. By using the central sub-images, says 3*3 or 4*4, we generate an image with deeper depth of field. We can also do some pre-process, to sample a larger gird on lens. When rendering, we just use the central grid so that after rendering we can also decrease the depth of field by increase the number of sub-images to generate image.
    </p>

    <p>
        Note that by extending depth of field in this way, we will get image with higher SNR because it wastes light of full aperture. To handle this problem, we could refocus each pixel to form the final image. This method will use all the information we received. However, since the accuracy of our depth detection algorithm is not high, we did not implement this approach.
    </p>

    <h4>4. Camera position adjustment</h4>

    <p>
        We use the same idea as depth of field adjustment in this part. In order to generate image with different camera positions, we use different sub-images. For example, if we want a image that taken by camera that is a little left than before, we can pick a sub-grid in the left of the full sub-images gird.
    </p>

    <p>
        The problem is that, since we use CPU to render the image, and need to do a lot pre-process, if we render a lager sub-image grid, it will take really long time and consume too much memory. While with a relatively small sample grid, the difference between each camera position is not so obvious. An approach to solve this problem is to use GPU to render instead CPU. However, we need to rewrite all the skeleton if we choose to render by GPU. Since we have limited time, we do not decided to use this method.
    </p>

    <h3 align="middle">Result</h3>

    <h4>1.The effect of refocusing.</h4>
    
    <h5>(1)CBdragon</h5>

    <div align="center">
    	<tr>
    		<br>
            	<img src="images/Refocus_dragon1.png" width="600px"/>
            	<figcaption align="middle">Focus at the head of the dragon.</figcaption>
        	</br>
        	<br>
            	<img src="images/Refocus_dragon2.png" width="600px"/>
            	<figcaption align="middle">Foucs at the body of the dragon.</figcaption>
        	</br>
        	<br>
            	<img src="images/Refocus_dragon.png" width="600px"/>
            	<figcaption align="middle">Focus at the front of the dragon, blurry.</figcaption>
        	</br>
        </tr>
    </div>

     <h5>(2)Rabbit</h5>

    <div align="center">
    	<tr>
    		<br>
            	<img src="images/Render_rabbit.png" width="600px"/>
            	<figcaption align="middle">Focus at rabbit.</figcaption>
        	</br>
        	<br>
            	<img src="images/Refocus_rabbit.png" width="600px"/>
            	<figcaption align="middle">Focus at the front of the rabbit, blurry.</figcaption>
        	</br>
        </tr>
    </div>




    <h4>2.Changing depth of field.</h4>

    <h5>(1)CBdragon</h5>

    <div align="center">
    	<tr>
    		<br>
            	<img src="images/Render_rabbit.png" width="600px"/>
            	<figcaption align="middle">Normal Aperture.</figcaption>
        	</br>
        	<br>
            	<img src="images/Aperture_increase_rabbit.png" width="600px"/>
            	<figcaption align="middle">Increase Aperture, small DOF, less noise.</figcaption>
        	</br>
        	<br>
            	<img src="images/Aperture_decrease_rabbit.png" width="600px"/>
            	<figcaption align="middle">Decrease Aperture, large DOF, more noise.</figcaption>
        	</br>
        </tr>
    </div>

    <h4>3.Camera postion adjustment.</h4>

    <h5>(1)CBdragon</h5>

    <div align="center">
    	<tr>
        	<br>
            	<img src="images/up.png" width="600px"/>
            	<figcaption align="middle">Shift lens up.</figcaption>
        	</br>
        	<br>
            	<img src="images/down.png" width="600px"/>
            	<figcaption align="middle">Shift lens down.</figcaption>
        	</br>
        	<br>
            	<img src="images/left.png" width="600px"/>
            	<figcaption align="middle">Shift lens left.</figcaption>
        	</br>
        	<br>
            	<img src="images/right.png" width="600px"/>
            	<figcaption align="middle">Shift lens right.</figcaption>
        	</br>
        </tr>
    </div>
    <!-- Resources -->
    <h3 align="middle">Reference</h3>
    <p>
        <a href="https://graphics.stanford.edu/papers/lfcamera/lfcamera-150dpi.pdf">https://graphics.stanford.edu/papers/lfcamera/lfcamera-150dpi.pdf</a>
    </p>
    <p>
        <a href="https://en.wikipedia.org/wiki/Light_field">https://en.wikipedia.org/wiki/Light_field</a>
    </p>
    <p>
        <a href="https://graphics.stanford.edu/papers/light/light-lores-corrected.pdf">https://graphics.stanford.edu/papers/light/light-lores-corrected.pdf</a>
    </p>

    <!-- Contributions -->
    <h3 align="middle">Contributions</h3>

    <h4>Renjie Shao</h4>
    <ul>
        <li>Implement light field information collection.</li>
        <li>Implement refocusing.</li>
        <li>Writing webpages.</li>
    </ul>
    
    <h4>Nan Wei</h4>
    <ul>
        <li>Implement depth of field adjustment.</li>
        <li>Implement position of camera adjustment.</li>
        <li>Writing webpages.</li>
    </ul>

    <h4>Jiyuan Lu</h4>

    <ul>
        <li>Initial the project.</li>
        <li>Present the project at the poster session.</li>
        <li>Make the milestone and presentation slides.</li>
        <li>Record the milestone and final project videos.</li>
    </ul>


<h3 align="Middle">Slides</h3>

<div align="middle">
    <iframe src="https://docs.google.com/presentation/d/e/2PACX-1vSEw8M2SAOPXlpG1ORFzmtNofZ96wXyOG3i_gZfB3mmE7hDEy4M4R-IW4D0irCharDErp_pG6PnLku6/embed?start=false&loop=false&delayms=3000" frameborder="0" width="1058" height="823" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>
</div>

<h3 align="middle">Video</h3>

<div align="middle">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/g7M8S4g8IAc" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
</div>

</div>
