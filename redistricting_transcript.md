# Redistricting Algorithm — Session Transcript & Outline

---

## Full Transcript

---

### User (Prompt 1)

Hey. So, um, I dictated this before, and it forgot it. So... and it literally took ten minutes. And you know what? I haven't had yet today lunch. So thank you, Claude, for wasting ten minutes of my time when I've been sort of on my computer for the past six hours working on this, but we'll get over it. So I'm gonna try to start from the beginning. So I like... this is... my first recording was so much better, by the way, Claude. So I hope you're happy that you just ruined my life. But I like doing spatial clustering and stuff like that, and I like, um, I've been having idea of, like, okay, redistricting algorithm. Fun. That's sort of similar to what I have been doing previously and, like, personal projects, whatever. So a few months ago, I heard about this algorithm called Metis, m e t i s, and what it is is, um, a clustering algorithm for, um, graphs, so node edge things. It is, shall we say, adept at creating min cuts of edges and, uh, partitioning graphs based on, um, some constraints. The constraint that is relevant for me here is, yeah, min cut and then ensuring that, um, the graph, um, what is it called, partitions are... have equal node rate. So that could just be number of nodes, or it could be, like, each node has a property and we're equalizing that. Um, obviously, for the problem of redistricting, node weight can be... it's still popularly shown. So, basically, we are trying to make min cuts that, um, redistrict places into equal population. Great. So why do this? Because it's something that the government does and something that's interesting is one that people talk about. Also, I... so what are, like, motivations? So there are basically, like, certain proper is that people like in any district in any district in the legislative district whatever. And there are some in particular that I personally like, um, that sometimes I think people don't always like or think about. I think usually people do, like, and think about them. I can't believe I had more oatmeal here, and I didn't know about it. Okay. That will be for tomorrow or tonight, though. Anyway, so compactness, obviously, the district to compete... be compact. You don't want, like, a gerrymandering, and you don't want, like, let's say, sometimes people don't like being disturbed by dividing it into... across, like, lines of latitude into... we'll divide New York into twenty six different districts. Okay. So that's annoying because then, like, Buffalo is going to be in the same district as somewhere in the border, like Massachusetts. It's just, like, a mess, and it's silly. Or, like, the smallest cut line where you just sort of, like, draw a straight line through the whole state until... and you keep doing that until you get what you want. It's like, okay. Well... and Noninsaneland, those don't sound like good ideas. Um, so the other thing... I'm gonna come back to this. It's, like, cohesiveness of, uh, communities. So, like, if you are in a city with seven hundred thousand people, probably you want them all to be in the same district because each congressional district in the US has close to seven hundred thousand people, a little more. So that's, uh, another thing. And one way to achieve this is to decide that you will put the district lines in sort of outlining areas where not a lot of people are living, so you're not affecting a lot. You're, like, not dividing up communities, um, as much. Like, for example, if you have a congressional district in New Jersey and it just goes, like, sort of straight through the middle of Newark, It's like, well, that's sort of not what you want because my work is sort of like a cohesive entity or, like... I mean, New York has several, so you have to have them go through, like, the middle of New York. Um, but, like, if you if you sort of had, like, a smaller city, you wouldn't want, like, two or even three lines. You just go through the middle of it. Maybe if there's a river, and that's actually what happens in that other mallet, which I think is a good thing. But, anyway, I'm getting beside the point. So a few weeks ago, I sort of just, like, experimenting with the Mentos algorithm or whatever, and I was like, okay. Well, I... it's cool that the... this does mean cuts, um, but, like, the parameters are not to be stone. Like, how do I decide? Do I use census tracts? Do I use census blocks? Do I use census plunk groups for deciding... and I see in my last recording, I was talking about how you can decide how places are even adjacent to each other, so I guess that is relevant. Like, you can decide, okay, census tracts. Census tracts are sort of, like, easy to talk about. You can say that they... so, like, you can arrange them into a graph. So you have, like, New York State has, like, maybe around three thousand census tracts. So you can you can basically say, okay. So if nice... if we're creating a graph out of census tracts, every node is going to be a centroid of a tract, and every edge is going to be between two tracts that border each other. So you could do it that way, or you could do it a different way. there are there are a few different ways too. Like, I I could do Delonate triangulation where you decide who is ordering what based on sort of a geometric optimization. Um, so it will not always be the case that, like, two census tracts that order either will, like, share an edge when you do the Delonate triangulation, but it sort of makes, like, an efficient mesh where I think, like, the circumcircle of any three points on a triangle. It creates a whole bunch of triangles that end up being the edges of the graph. The circle that circumscribes that triangle will not ever include any other nodes. So that's sort of an efficient way of just creating nodes or creating edges for nodes where you, like, don't know, like, there's no really objective, like, way other than, like, distance to determine who borders what or what is near what.

---

### Claude (Response 1)

So I like doing spatial clustering, and I've been toying with an idea: a redistricting algorithm. It's similar to personal projects I've worked on before.

A few months ago, I heard about an algorithm called METIS — M-E-T-I-S. It's a clustering algorithm for graphs (node-edge structures). It's particularly good at creating minimum cuts of edges and partitioning graphs based on constraints. The constraint most relevant here is the min-cut combined with ensuring that graph partitions have equal node weight. Node weight can just be the number of nodes, or each node can carry a property that gets equalized across partitions. For redistricting, that property is population — we're trying to make minimum cuts that divide places into districts of roughly equal population.

**Why do this?** Redistricting is something governments actually do, it's something people care about, and there are interesting properties that any good legislative district should have.

**Compactness** is the obvious one. You don't want gerrymandering, and you don't want absurd configurations — like dividing New York into 26 horizontal slices so that Buffalo ends up in the same district as somewhere near the Massachusetts border. Or just drawing straight lines through a state until you get the number of districts you want. Those approaches produce nonsensical maps.

**Community cohesion** is another property I care about. If you have a city of 700,000 people, you probably want them all in the same congressional district — since each US congressional district is roughly that size anyway. One way to achieve this is to draw district lines through low-population areas, so you're not cutting through communities. If a congressional district line goes straight through the middle of Newark, that's not ideal. Newark is a cohesive entity. For a smaller city, you'd ideally avoid running two or three lines through the center of it. (Sometimes a river provides a natural boundary, and that can actually work well.)

---

A few weeks ago I started experimenting with METIS, and pretty quickly ran into some open questions. The algorithm handles the min-cuts, but the parameters aren't set in stone. What geographic units do you use — census tracts, census blocks, block groups? And how do you decide which units are even *adjacent* to each other?

Census tracts are a reasonable starting point. New York State has roughly 3,000 of them. You can construct a graph where every node is the centroid of a tract, and every edge connects two tracts that share a border. That's one approach.

Another option is **Delaunay triangulation** — a geometric method for deciding adjacency. It won't always match actual shared borders, but it produces an efficient mesh. The defining property is that for any triangle formed by three nodes, the circumscribed circle around that triangle contains no other nodes. It's a clean way to define "who is near what" when there's no obvious objective criterion beyond distance.

---

### User (Prompt 2)

So, Mike, this is great. Um, I can literally just run lettuce on, uh, adjacency matrix of census tracts and get pretty good districts. Um, I spent a few more days experimenting with, like, other things. Like, I don't know. What was it? Well, now I've distracted myself late, like, actually thinking of, like, what's happening, and I'm gonna try to undistract myself. So, yeah, I was like, okay. Well, I'm not gonna use such distractions. Those are not, like, really cool. Probably, like, two thousand to eight thousand is a big difference. I'm going to just buy myself a group together and relax and try to have an equal population. I was doing that with hexagon meshes. I was doing that with the meshes algorithm. I was doing that with... I think it was a major way to those who might get satisfactory results because they kept, like, getting disconnected components and things. Like, I I I would bring the box together, and I would be, like, not... they shouldn't have knocked you out. What I ended up doing was literally just going back to the center. Like, I spent, like, maybe two more afternoons on that. I was like, And then I was like, well, what if I just go back to census tracts? And that, you know, back to basics. And that was a really good outcome. And I'm like, okay. But some census checks are really still too large. And really, the issue with using census checks is that they include water. And if you are, I mean, someone from Long Island, like, the water... just because something is as the crow flies at distance, that doesn't mean that they're actually connected. Under seats that are like that, Massachusetts to an extent, Rhode Island and Maryland to a large extent, Virginia, sort of Michigan. These are states that have, like, a lot of islands and just areas that are rated by water that you don't want to just assume are close to each other just because they might be, like, two miles away. Maybe two miles away and very far from any populated place or bridge or ferry service. So that was something to keep in mind. And if you go to the blog level. Mhmm. Water blocks are always separate from land blocks, which is not... for the census block level, often water is separate from, uh, land and at the census tract level, it's pretty rare that water is separate from land. I think only for, like, the largest libraries of water. So, yeah, I don't want I don't want it to just be assumed that, like, Glen Cove borders, like, City Island and the Bronx. Yeah. So I did end up using... so I'm just trying to as a base, but I have two sort of constraints. Two constraints as far as, like, the size of the nodes that I'm putting in. So one constraint is that -- Mhmm. -- if taking all the box of a... all the populated mocks of a census tract makes the tract discontinuous, then considered as if it's two or more separate tracts, um, that don't touch each other. So that's one thing. And another thing is that if a centist tract just has, like, a really high population, split it into smaller tracts that have less than two times the median tract population. Another thing is that if you have a bunch of disconnected components, um, merge them there's a whole lot of logic that I had to sort of decide how to merge disconnect components that I won't get into because it's annoying. And I might want to refine it. Like, for... I I think I've... now we've done a pretty good job of, like, saying, connect islands to the mainland that don't do... don't sort of assume that they're, like, connected in the same way that anyone else would be. Like, maybe cut the island off and put it in a different district if it makes sense. And then there's a whole issue of, like, Long Island being connected to, like, Westchester for some reason. Or, like, I even had, like, the area of New York close to close to where Quebec begins, like, on the western side of, like, the Canada border of New York being directly connected to Plattsburgh. And I'm like, that shouldn't happen. All of these... so I think I will probably... the next step will be to, like, have inputs that are like... when you decide that tracts border each other and when you draw lines between disconnected components, do not draw a line through this line. Like, do not draw a line through the lung sound because we are not we are not having, like, Huntington me in the same congressional district as Poughkeepsie. It's just not gonna happen. So... and, really, like, the line should only rarely go the outside of the border if the sand was the border of the state is, like, so large, which I think is not Something that happens on this floor today. Oh, it does happen for Nantucket. Nantucket is like... anyway, like, there there needs to be a few more moments, but I haven't done that yet. So we have basically census tract shaped errors. Like, some of the population now we have closer and... all closer in population. You do make cuts, and they're fairly close to the census. Um, what is it called? Population... the congressional district population. They are reasonable. I think if you're using pure straights, they are, like, very reasonable in shape. When I'm, like, spinning them up, sometimes they get a little wonky. Oh, and also because, like, the senses tracks will... there's no, like, dealing with discontent components because of the water. So that's, like, easier. it does get a little wonky at the census tract. Split it in a water issue. Some cities do get split in half because they're sort of popular stuff that they cannot be all in one. And so at that point, it's like... yeah. So the whole, like, Ontario, west of Rochester area, downtown Tampa, Buffalo is one district. And then the whole other half of Buffalo, it dumps, like, everywhere sort of south of it in New York state as another district. And it's like, I could... like, doing a mooncut cutting out, like, a large amount of silver in here in county would... that would not be actually conducive to a rain cut because of a a large number of somewhat populous sets of strands there or a large number of sets of strands in there. So as you cut right through Buffalo, you are sort of cutting through for a short time through, like, Buffalo, but then you're, like, out in the rural area. So it's basically, like, Buffalo North, Buffalo South.

---

### Claude (Response 2)

So the good news: you can literally just run METIS on an adjacency matrix of census tracts and get pretty reasonable districts. But I spent a few more days experimenting with alternatives — hexagon meshes, other binning approaches — and kept running into the same problem: disconnected components. I'd group things together and end up with clusters that weren't geographically contiguous. After a couple of afternoons on that, I went back to basics. Census tracts. And honestly, that gave the best results.

That said, census tracts come with their own issues.

**The water problem**

Census tracts typically include water area, which means two tracts that are only a couple miles apart as the crow flies might be separated by a bay, a sound, or a river with no bridge or ferry nearby. States where this really matters: Maryland, Rhode Island, Virginia, Michigan — and New York, where you really don't want to assume that Glen Cove on Long Island borders City Island in the Bronx just because of proximity. At the census block level, water and land are usually kept separate. At the tract level, that separation only happens for the largest bodies of water, so it's something you have to handle explicitly.

**Preprocessing the tracts**

To address this, I added two constraints on the nodes going into the graph:

1. **Strip water, check connectivity.** If removing a tract's water-only blocks makes it geographically discontinuous, treat it as two or more separate tracts that don't touch each other.
2. **Split oversized tracts.** If a tract's population is unusually high, split it into sub-tracts each with less than twice the median tract population.

There's also a fair amount of logic for handling disconnected components — merging them in sensible ways. I won't get into all of it, but the goal is to connect islands to the mainland without pretending they're just like any other neighboring tract. Sometimes it makes more sense to cut an island off and put it in a different district entirely.

Some edge cases still need work. I ended up with Long Island connecting to Westchester in ways it shouldn't, and at one point the area of upstate New York near the Quebec border was getting linked directly to Plattsburgh in a way that made no geographic sense. The next step is probably adding explicit barrier inputs — lines the algorithm is not allowed to draw edges across. No line through Long Island Sound. No edge that puts Huntington and Poughkeepsie in the same district. The border of a state should almost never be crossed, with maybe a rare exception like Nantucket. Those guardrails aren't built yet.

**Results so far**

With the preprocessed census tracts, population is now much more consistent across nodes, and the METIS cuts produce districts that are fairly close to the target congressional population. For purely inland states, the shapes are actually quite reasonable.

Where it gets a little wonky is in areas affected by the water and split-tract logic. Some cities do end up divided — when a city is populous enough that it can't all fit in one congressional district, the algorithm cuts straight through it. Buffalo is a good example: one district ends up being roughly the northern half, another gets the southern half plus a large swath of rural western New York. That happens because cutting across Buffalo — even though you're briefly passing through a dense area — quickly gets you into lower-population rural tracts, so it's actually an efficient min-cut. It's not ideal, but it's not unreasonable either. "Buffalo North, Buffalo South" is at least a coherent outcome.

---

### User (Prompt 3)

So, Mike, this is great. Um, I can literally just run lettuce on, uh, adjacency matrix of census tracts and get pretty good districts. Um, I spent a few more days experimenting with, like, other things. Like, I don't know. What was it? Well, now I've distracted myself late, like, actually thinking of, like, what's happening, and I'm gonna try to undistract myself. So, yeah, I was like, okay. Well, I'm not gonna use such distractions. Those are not, like, really cool. Probably, like, two thousand to eight thousand is a big difference. I'm going to just buy myself a group together and relax and try to have an equal population. I was doing that with hexagon meshes. I was doing that with the meshes algorithm. I was doing that with... I think it was a major way to those who might get satisfactory results because they kept, like, getting disconnected components and things. Like, I I I would bring the box together, and I would be, like, not... they shouldn't have knocked you out. What I ended up doing was literally just going back to the center. Like, I spent, like, maybe two more afternoons on that. I was like, And then I was like, well, what if I just go back to census tracts? And that, you know, back to basics. And that was a really good outcome. And I'm like, okay. But some census checks are really still too large. And really, the issue with using census checks is that they include water. And if you are, I mean, someone from Long Island, like, the water... just because something is as the crow flies at distance, that doesn't mean that they're actually connected. Under seats that are like that, Massachusetts to an extent, Rhode Island and Maryland to a large extent, Virginia, sort of Michigan. These are states that have, like, a lot of islands and just areas that are rated by water that you don't want to just assume are close to each other just because they might be, like, two miles away. Maybe two miles away and very far from any populated place or bridge or ferry service. So that was something to keep in mind. And if you go to the blog level. Mhmm. Water blocks are always separate from land blocks, which is not... for the census block level, often water is separate from, uh, land and at the census tract level, it's pretty rare that water is separate from land. I think only for, like, the largest libraries of water. So, yeah, I don't want I don't want it to just be assumed that, like, Glen Cove borders, like, City Island and the Bronx. Yeah. So I did end up using... so I'm just trying to as a base, but I have two sort of constraints. Two constraints as far as, like, the size of the nodes that I'm putting in. So one constraint is that -- Mhmm. -- if taking all the box of a... all the populated mocks of a census tract makes the tract discontinuous, then considered as if it's two or more separate tracts, um, that don't touch each other. So that's one thing. And another thing is that if a centist tract just has, like, a really high population, split it into smaller tracts that have less than two times the median tract population. Another thing is that if you have a bunch of disconnected components, um, merge them there's a whole lot of logic that I had to sort of decide how to merge disconnect components that I won't get into because it's annoying. And I might want to refine it. Like, for... I I think I've... now we've done a pretty good job of, like, saying, connect islands to the mainland that don't do... don't sort of assume that they're, like, connected in the same way that anyone else would be. Like, maybe cut the island off and put it in a different district if it makes sense. And then there's a whole issue of, like, Long Island being connected to, like, Westchester for some reason. Or, like, I even had, like, the area of New York close to close to where Quebec begins, like, on the western side of, like, the Canada border of New York being directly connected to Plattsburgh. And I'm like, that shouldn't happen. All of these... so I think I will probably... the next step will be to, like, have inputs that are like... when you decide that tracts border each other and when you draw lines between disconnected components, do not draw a line through this line. Like, do not draw a line through the lung sound because we are not we are not having, like, Huntington me in the same congressional district as Poughkeepsie. It's just not gonna happen. So... and, really, like, the line should only rarely go the outside of the border if the sand was the border of the state is, like, so large, which I think is not Something that happens on this floor today. Oh, it does happen for Nantucket. Nantucket is like... anyway, like, there there needs to be a few more moments, but I haven't done that yet. So we have basically census tract shaped errors. Like, some of the population now we have closer and... all closer in population. You do make cuts, and they're fairly close to the census. Um, what is it called? Population... the congressional district population. They are reasonable. I think if you're using pure straights, they are, like, very reasonable in shape. When I'm, like, spinning them up, sometimes they get a little wonky. Oh, and also because, like, the senses tracks will... there's no, like, dealing with discontent components because of the water. So that's, like, easier. it does get a little wonky at the census tract. Split it in a water issue. Some cities do get split in half because they're sort of popular stuff that they cannot be all in one. And so at that point, it's like... yeah. So the whole, like, Ontario, west of Rochester area, downtown Tampa, Buffalo is one district. And then the whole other half of Buffalo, it dumps, like, everywhere sort of south of it in New York state as another district. And it's like, I could... like, doing a mooncut cutting out, like, a large amount of silver in here in county would... that would not be actually conducive to a rain cut because of a a large number of somewhat populous sets of strands there or a large number of sets of strands in there. So as you cut right through Buffalo, you are sort of cutting through for a short time through, like, Buffalo, but then you're, like, out in the rural area. So it's basically, like, Buffalo North, Buffalo South.

---

### Claude (Response 3)

Here's the cleaned-up continuation for your blog:

---

So there are several parameters you can tune when setting up this kind of redistricting graph. I've already covered two: what the nodes are, and how you define adjacency (i.e., where the edges go). Let me go deeper on both, and then get into edge weights.

**On nodes:** Census tracts are a natural starting point. As background — the Census Bureau tries to keep tracts between 2,000 and 8,000 people. Block groups sit below that, targeting roughly 300 to 800 people. And then census blocks are the finest level — literally a city block, or any discrete parcel of land that's easy to identify. There's no population target at the block level; it could be a dense urban block or a sprawling area of cul-de-sacs that don't even share a street but get treated as one unit.

I was also experimenting with an alternative: binning census blocks into hexagons, and if a hexagon's population got too large, splitting it into smaller hexagons. Node weight in that case is just the sum of population within each bin.

I'm running this on a personal computer, so scale matters. I *can* run METIS on every census block, but it gets slow — and what really takes time is the Delaunay triangulation, because census tracts often share only a small sliver of a border, and accounting for that gets expensive. My practical workaround: only give two tracts an edge if they physically border each other *and* the Delaunay triangulation also gives them an edge. That combination filters out a lot of noise.

**On edge weights:** This is where things got interesting.

METIS does minimum cuts, so a high edge weight means "expensive to cut here" — i.e., you want to keep those nodes together. A low edge weight means the algorithm will happily cut that edge.

My first instinct was to use **1/distance** as the edge weight. Far-apart tracts get a low weight, so they're easy to cut. Close tracts get a high weight, so they stay together. Intuitive.

I also tried weighting by population — something like **distance / √population**, which is roughly proportional to the average distance between two people across the edge. High population nodes that are close together get an expensive edge, which discourages cutting through dense areas.

Neither of those gave me the results I wanted. They tended to just cluster adjacent tracts that happened to border each other.

What actually worked better? **Uniform edge weights** — every edge gets a weight of 1. I did not expect that to work at all. But here's why it makes sense:

Census tracts are designed to be roughly equal in population. That means low-population tracts tend to be geographically large (think rural areas), and high-population tracts tend to be small (dense urban areas). When you do a min-cut with uniform edge weights, cutting a single edge between two large rural tracts removes a big chunk of geographic area from a subgraph in one move — you've made a lot of progress toward a valid partition with just one cut.

In other words, the geographic size of low-population tracts already encodes the population information implicitly. I don't need to manually weight by √population over distance — that effect is already baked into the graph structure. The uniform weight approach naturally steers cuts toward sparse, large-area regions, which is exactly what good redistricting should do: draw lines where fewer people live.

---

### User (Prompt 4)

You know, if other results are, like, Manhattan usually ends up in parts of three districts, so, like, south, middle, north, western part of the Bronx is sometimes... yeah. Usually, there's western part of the Bronx district and eastern part of the Bronx district. Um, there's the districts, like, Westchester and then, like, maybe Orange and Rockland and then, like, the area. Like, the new districts are pretty much like that, and that's about We have Pacific County, Easton Park in the district. You have Western, Suffolk, and two, Nassau District. You have a North Shore, Nassau District, and a South Shore, Nassau District. Yeah. So... oh, and you have the statin in the district that includes a reasonable close area of Brooklyn. So I think that the district has gone up. There have been a lot of times where I'm like, okay. So paper cut is not part of the North Shore of Massachusetts actually. So why are they in this industry? Are -- I think I'm just tweaking for it to not be like that. I will not do that. I just wanna make a freaking blog about it. I don't wanna spend this whole week on that. I wanna do the blog and put it out and then maybe do a revision later. Thank you. Thank you, and good night.

---

### Claude (Response 4)

**Results for New York**

The New York districts the algorithm produces are actually pretty recognizable. Manhattan typically ends up split into three districts — roughly south, middle, and north. The Bronx usually divides into a western district and an eastern district. Westchester gets its own district, often paired with something like Orange and Rockland counties. Long Island breaks down into Western Suffolk, two Nassau districts (one north shore, one south shore), and an Eastern Suffolk/East End district. Staten Island ends up in a district that includes a reasonable adjacent slice of Brooklyn.

Overall, the shapes are pretty sensible. There are occasional oddities — a town ending up on the wrong side of a boundary, where you look at it and think "why is this place grouped with that?" — but rather than keep tweaking forever, the plan is to write this up, put it out, and revisit later.

---

## Outline of What Was Said

### 1. Introduction & Motivation
- Speaker is interested in spatial clustering and redistricting algorithms
- Discovered METIS a few months ago: a graph partitioning algorithm
- METIS specializes in minimum edge cuts with constraints, particularly equal node weight
- For redistricting, node weight = population; goal is equal-population districts
- Redistricting is something real governments do and people care about

### 2. Desirable Properties of Districts
- **Compactness:** No gerrymandering; no absurd horizontal slicing of states; no simple straight-line cuts
- **Community cohesion:** Cities should stay together where possible; district lines should run through low-population areas; rivers can serve as natural boundaries

### 3. Graph Construction Parameters
- **Node choice:** Census tracts (2,000–8,000 people), block groups (300–800), or census blocks (no population target, just discrete land parcels)
- Also experimented with hexagon bins — grouping census blocks into hexagons and splitting high-population ones; node weight = sum of population
- **Adjacency / edge definition:**
  - Shared-border approach: edge between any two tracts that physically border each other
  - Delaunay triangulation: geometric method; circumscribed circle of any triangle contains no other nodes; efficient mesh for defining proximity without obvious objective criterion
  - Practical hybrid: only draw an edge if tracts both physically border AND triangulation agrees

### 4. Edge Weights — What Worked and What Didn't
- **1/distance:** Far tracts easy to cut, close tracts expensive. Intuitive but didn't give good results.
- **distance / √population:** Approximates average distance between people across an edge. Also unsatisfying — tended to cluster only immediately adjacent tracts.
- **Uniform weight (1):** Counterintuitive but worked best. Reason: low-population tracts are geographically large, so cutting one edge between rural tracts covers a lot of geographic ground. Population information is already implicit in tract size — no manual weighting needed.

### 5. Node Preprocessing — Handling Water and Oversized Tracts
- Census tracts include water, which creates false adjacency (e.g., Glen Cove appearing to border City Island)
- States most affected: Maryland, Rhode Island, Virginia, Michigan, New York
- At block level, water is usually separated from land; at tract level, only the largest water bodies are separated
- **Preprocessing rules applied:**
  1. Strip water blocks; if remaining land is discontinuous, split into separate tracts
  2. Split high-population tracts into sub-tracts under 2x median tract population
  3. Merge disconnected components using custom logic (not fully detailed)
- **Known remaining issues:**
  - Long Island connecting to Westchester inappropriately
  - Upstate NY near Quebec border connecting directly to Plattsburgh
- **Planned fix:** Explicit barrier lines — edges cannot cross Long Island Sound, state borders, etc.

### 6. Results
- **General:** Districts are fairly close to congressional population targets; shapes are reasonable for inland areas
- **City splitting:** Dense cities sometimes get divided (e.g., Buffalo → "Buffalo North" and "Buffalo South") because cutting through a dense city briefly is still an efficient min-cut when it quickly reaches rural tracts
- **New York City specifically:**
  - Manhattan: three districts (south, middle, north)
  - Bronx: western and eastern districts
  - Westchester: paired with Orange and Rockland counties
  - Long Island: Western Suffolk, two Nassau (north shore / south shore), Eastern Suffolk
  - Staten Island: paired with adjacent part of Brooklyn
- Some oddities remain (unexpected town groupings) but the plan is to publish now and revise later
